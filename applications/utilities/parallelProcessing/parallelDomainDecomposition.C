/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "parallelDomainDecomposition.H"
#include "cpuTime.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "fvCFD.H"

// change the constructor
Foam::parallelDomainDecomposition::parallelDomainDecomposition(
    const IOobject &io,
    label procNo)
    : fvMesh(io),
      facesInstancePointsPtr_(
          pointsInstance() != facesInstance()
              ? new pointIOField(
                    IOobject(
                        "points",
                        facesInstance(),
                        polyMesh::meshSubDir,
                        *this,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false))
              : nullptr),
      procNo_(
          procNo),
      cellToProc_(nCells()),
      procPointAddressing_(8),
      procFaceAddressing_(8),
      procCellAddressing_(8),
      procPatchSize_(8),
      procPatchStartIndex_(8),
      procNeighbourProcessors_(8),
      procProcessorPatchSize_(8),
      procProcessorPatchStartIndex_(8),
      procProcessorPatchSubPatchIDs_(8),
      procProcessorPatchSubPatchStarts_(8)
{
    // TODO: change constructor initializer list
    Pout << "parallelDomainDecomposition constructor is called for process " << procNo_ << endl;
}

void Foam::parallelDomainDecomposition::decomposeMesh()
{
    Pout << "decomposeMesh method is called" << endl;
    distributeCells();
}

bool Foam::parallelDomainDecomposition::writeDecomposition()
{
    Pout << "writeDecomposition is called" << endl;
    Pout << "\nCalculating original mesh data for process " << Pstream::myProcNo() << endl;

    return true;
}

void Foam::parallelDomainDecomposition::distributeCells()
{
    Info<< "\nCalculating distribution of cells" << endl;

    cpuTime decompositionTime;
    const decompositionModel& method = decompositionModel::New
    (
        *this
    );

    word weightName;
    scalarField cellWeights;

    if (method.readIfPresent("weightField", weightName))
    {
        volScalarField weights
        (
            IOobject
            (
                weightName,
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            *this
        );
        cellWeights = weights.primitiveField();
    }

    cellToProc_ = method.decomposer().decompose(*this, cellWeights);

    // Pout << cellToProc_ << endl;
    // autoPtr<OFstream> outputFilePtr;
    // fileName outputDir = "/home/dogukan/Documents/hydrology/tutorials/permaFoam/demoCase/debug";
    // outputFilePtr.reset(new OFstream(outputDir/"cellToProc_"));
    // outputFilePtr() << cellToProc_ << endl;

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}

// ************************************************************************* //
