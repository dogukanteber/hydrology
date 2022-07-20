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

    // get number of processes from the dict
    scalar nProc(readScalar(decompDict.lookup("numberOfSubdomains")));

    Info<< "Create undecomposed database"<< nl << endl;
    Time baseRunTime
    (
        runTime.controlDict(),
        runTime.rootPath(),
        runTime.globalCaseName(),
        runTime.system(),
        runTime.constant(),
        false                   // enableFunctionObjects
    );

    Pout << "Create mesh for process " << Pstream::myProcNo() << endl;
    parallelDomainDecomposition mesh
    (
        IOobject
        (
            "region0",  // change this later!!
            baseRunTime.timeName(),
            baseRunTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Pstream::myProcNo()
    );

    mesh.decomposeMesh();
    mesh.writeDecomposition();

    return 0;
}

// ************************************************************************* //
