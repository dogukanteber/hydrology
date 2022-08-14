/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    parallelDecompose

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case in parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "error.H"
#include "parallelDomainDecomposition.H"
#include "fvFieldDecomposer.H"
#include "IOobjectList.H"
#include "readFields.H"
#include "pointFieldDecomposer.H"

/*
    Possible implementations:

        -> Add --quite option to not print output to the console (might be useful with many processes)
        -> Return an error if the mesh is already decomposed
        -> Return an error if blockMesh is not called
*/


// Read proc addressing at specific instance.
// Uses polyMesh/fvMesh meshSubDir by default
autoPtr<labelIOList> procAddressing
(
    const fvMesh& procMesh,
    const word& name,
    const word& instance,
    const word& local = fvMesh::meshSubDir
)
{
    return autoPtr<labelIOList>::New
    (
        IOobject
        (
            name,
            instance,
            local,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // do not register
        )
    );
}

// Return cached or read proc addressing from facesInstance
const labelIOList& procAddressing
(
    const fvMesh& procMesh,
    const word& name,
    PtrList<labelIOList>& procAddressingList,
    const label proci = 0
)
{
    if (!procAddressingList.set(proci))
    {
        procAddressingList.set
        (
            proci,
            procAddressing(procMesh, name, procMesh.facesInstance())
        );
    }
    return procAddressingList[proci];
}


int main(int argc, char* argv[]) {

    argList::addNote
    (
        "Decompose a mesh and fields of a case in parallel execution"
    );
    argList::noCheckProcessorDirectories();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Use specified file for decomposePar dictionary"
    );
    
    fileHandler().distributed(true);

    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << ": This utility can only be run parallel"
            << exit(FatalError);
    }

    #include "createTime.H"
    runTime.functionObjects().off();  // Extra safety?

    // get custom decomposeParDict location
    fileName decompDictFile(args.getOrDefault<fileName>("decomposeParDict", ""));
    if (!decompDictFile.empty() && !decompDictFile.isAbsolute())
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }
    else {
        decompDictFile = "decomposeParDict";
    }

    // create a dictionary object
    IOdictionary decompDict
    (
        IOobject
        (
            decompDictFile,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    // get decomposition method
    word method;
    decompDict.lookup("method") >> method;

    // if method is not simple, print error and exit the program
    if (method != "simple") {
        FatalErrorInFunction
                    << "Specified method type " << method
                    << " is not implemented"
                    << exit(FatalError);
    }

    Time baseRunTime
    (
        runTime.controlDict(),
        runTime.rootPath(),
        runTime.globalCaseName(),
        runTime.system(),
        runTime.constant(),
        false                   // enableFunctionObjects
    );

    instantList times = timeSelector::selectIfPresent(baseRunTime, args);

    parallelDomainDecomposition mesh
    (
        IOobject
        (
            parallelDomainDecomposition::polyMesh::defaultRegion,
            baseRunTime.timeName(),
            baseRunTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        decompDict,
        Pstream::myProcNo(),
        Pstream::nProcs()
    );

    // #include "PstreamReduceOps.H"
    // labelList data(8, 0);
    // data[Pstream::myProcNo() * 2] = Pstream::myProcNo();
    // data[Pstream::myProcNo() * 2 + 1] = Pstream::myProcNo();

    // reduce(data, sumOp<labelList>());

    // Info << data << endl;

    mesh.decomposeMesh();
    mesh.writeDecomposition();

    fileHandler().flush();

    // Decompose the field files

    // Cached processor meshes and maps. These are only preserved if
    // running with multiple times.

    // Since each process has its own memory, below lists are initialized
    // with only one element. Later in the code, you will see the relevant
    // data is access using Zero keyword. So, do not confuse.
    PtrList<labelIOList> faceProcAddressingList(One);
    PtrList<labelIOList> cellProcAddressingList(One);
    PtrList<labelIOList> boundaryProcAddressingList(One);
    PtrList<fvFieldDecomposer> fieldDecomposerList(One);
    PtrList<labelIOList> pointProcAddressingList(One);
    PtrList<pointFieldDecomposer> pointFieldDecomposerList
    (
        One
    );

    // Loop over all times
    forAll(times, timeI)
    {
        baseRunTime.setTime(times[timeI], timeI);

        Info<< "Time = " << baseRunTime.timeName() << endl;

        // Search for list of objects for this time
        IOobjectList objects(mesh, baseRunTime.timeName());

        // Construct the vol fields
        // ~~~~~~~~~~~~~~~~~~~~~~~~
        PtrList<volScalarField> volScalarFields;
        readFields(mesh, objects, volScalarFields, false);
        PtrList<volVectorField> volVectorFields;
        readFields(mesh, objects, volVectorFields, false);
        PtrList<volSphericalTensorField> volSphericalTensorFields;
        readFields(mesh, objects, volSphericalTensorFields, false);
        PtrList<volSymmTensorField> volSymmTensorFields;
        readFields(mesh, objects, volSymmTensorFields, false);
        PtrList<volTensorField> volTensorFields;
        readFields(mesh, objects, volTensorFields, false);


        // Construct the dimensioned fields
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
        readFields(mesh, objects, dimScalarFields);
        PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
        readFields(mesh, objects, dimVectorFields);
        PtrList<DimensionedField<sphericalTensor, volMesh>>
            dimSphericalTensorFields;
        readFields(mesh, objects, dimSphericalTensorFields);
        PtrList<DimensionedField<symmTensor, volMesh>>
            dimSymmTensorFields;
        readFields(mesh, objects, dimSymmTensorFields);
        PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;
        readFields(mesh, objects, dimTensorFields);


        // Construct the surface fields
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PtrList<surfaceScalarField> surfaceScalarFields;
        readFields(mesh, objects, surfaceScalarFields, false);
        PtrList<surfaceVectorField> surfaceVectorFields;
        readFields(mesh, objects, surfaceVectorFields, false);
        PtrList<surfaceSphericalTensorField>
            surfaceSphericalTensorFields;
        readFields(mesh, objects, surfaceSphericalTensorFields, false);
        PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
        readFields(mesh, objects, surfaceSymmTensorFields, false);
        PtrList<surfaceTensorField> surfaceTensorFields;
        readFields(mesh, objects, surfaceTensorFields, false);


        // Construct the point fields
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~
        const pointMesh& pMesh = pointMesh::New(mesh);

        PtrList<pointScalarField> pointScalarFields;
        readFields(pMesh, objects, pointScalarFields, false);
        PtrList<pointVectorField> pointVectorFields;
        readFields(pMesh, objects, pointVectorFields, false);
        PtrList<pointSphericalTensorField> pointSphericalTensorFields;
        readFields(pMesh, objects, pointSphericalTensorFields, false);
        PtrList<pointSymmTensorField> pointSymmTensorFields;
        readFields(pMesh, objects, pointSymmTensorFields, false);
        PtrList<pointTensorField> pointTensorFields;
        readFields(pMesh, objects, pointTensorFields, false);

        Pout << "Processor " << mesh.procNo() << ": field transfer" << endl;

        // open the database
        Time tempDb(
            Time::controlDictName,
            baseRunTime.rootPath(),
            baseRunTime.caseName() / ("processor" + Foam::name(mesh.procNo()))
        );

        Time& processorDb = tempDb;

        // read the mesh
        fvMesh tempMesh
        (
            IOobject
            (
                parallelDomainDecomposition::polyMesh::defaultRegion,
                processorDb.timeName(),
                processorDb
            )
        );
        const fvMesh& procMesh = tempMesh;

        const labelIOList& faceProcAddressing = procAddressing
        (
            procMesh,
            "faceProcAddressing",
            faceProcAddressingList
        );

        const labelIOList& cellProcAddressing = procAddressing
        (
            procMesh,
            "cellProcAddressing",
            cellProcAddressingList
        );

        const labelIOList& boundaryProcAddressing = procAddressing
        (
            procMesh,
            "boundaryProcAddressing",
            boundaryProcAddressingList
        );

        // FV fields: volume, surface, internal
        {
            if (!fieldDecomposerList.set(Zero))
            {
                fieldDecomposerList.set
                (
                    Zero,
                    new fvFieldDecomposer
                    (
                        mesh,
                        procMesh,
                        faceProcAddressing,
                        cellProcAddressing,
                        boundaryProcAddressing
                    )
                );
            }
            const fvFieldDecomposer& fieldDecomposer =
                fieldDecomposerList[Zero];

            // vol fields
            fieldDecomposer.decomposeFields(volScalarFields);
            fieldDecomposer.decomposeFields(volVectorFields);
            fieldDecomposer.decomposeFields
            (
                volSphericalTensorFields
            );
            fieldDecomposer.decomposeFields(volSymmTensorFields);
            fieldDecomposer.decomposeFields(volTensorFields);

            // surface fields
            fieldDecomposer.decomposeFields(surfaceScalarFields);
            fieldDecomposer.decomposeFields(surfaceVectorFields);
            fieldDecomposer.decomposeFields
            (
                surfaceSphericalTensorFields
            );
            fieldDecomposer.decomposeFields
            (
                surfaceSymmTensorFields
            );
            fieldDecomposer.decomposeFields(surfaceTensorFields);

            // internal fields
            fieldDecomposer.decomposeFields(dimScalarFields);
            fieldDecomposer.decomposeFields(dimVectorFields);
            fieldDecomposer.decomposeFields(dimSphericalTensorFields);
            fieldDecomposer.decomposeFields(dimSymmTensorFields);
            fieldDecomposer.decomposeFields(dimTensorFields);

            if (times.size() == 1)
            {
                // Clear cached decomposer
                fieldDecomposerList.set(Zero, nullptr);
            }
        }

        // Point fields
        if
        (
            pointScalarFields.size()
            || pointVectorFields.size()
            || pointSphericalTensorFields.size()
            || pointSymmTensorFields.size()
            || pointTensorFields.size()
        )
        {
            const labelIOList& pointProcAddressing = procAddressing
            (
                procMesh,
                "pointProcAddressing",
                pointProcAddressingList
            );

            const pointMesh& procPMesh = pointMesh::New(procMesh);

            if (!pointFieldDecomposerList.set(Zero))
            {
                pointFieldDecomposerList.set
                (
                    Zero,
                    new pointFieldDecomposer
                    (
                        pMesh,
                        procPMesh,
                        pointProcAddressing,
                        boundaryProcAddressing
                    )
                );
            }
            const pointFieldDecomposer& pointDecomposer =
                pointFieldDecomposerList[Zero];

            pointDecomposer.decomposeFields(pointScalarFields);
            pointDecomposer.decomposeFields(pointVectorFields);
            pointDecomposer.decomposeFields
            (
                pointSphericalTensorFields
            );
            pointDecomposer.decomposeFields(pointSymmTensorFields);
            pointDecomposer.decomposeFields(pointTensorFields);


            if (times.size() == 1)
            {
                pointProcAddressingList.set(Zero, nullptr);
                pointFieldDecomposerList.set(Zero, nullptr);
            }
        }
    }


    return 0;
}

// ************************************************************************* //
