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
          procNo)
{
    // TODO: change constructor initializer list
    Pout << "parallelDomainDecomposition constructor is called for process " << procNo_ << endl;
}

void Foam::parallelDomainDecomposition::decomposeMesh()
{
    Pout << "decomposeMesh method is called" << endl;
    distributeCells();

    Info<< "\nCalculating original mesh data" << endl;

    // set references to the original mesh
    const polyBoundaryMesh& patches = boundaryMesh();
    const faceList& fcs = faces();
    const labelList& owner = faceOwner();
    const labelList& neighbour = faceNeighbour();

    fileName path("/home/dogukan/Documents/hydrology/tutorials/permaFoam/demoCase/debugprocessor" + Foam::name(procNo_));
    mkDir(path);
    autoPtr<OFstream> outputFilePtr;

    Info<< "\nDistributing cells to processors" << endl;
    procCellAddressing_ = assignCellsToProc(cellToProc_);

    outputFilePtr.reset(new OFstream(path/"procCellAddressing_"));
    outputFilePtr() << procCellAddressing_ << endl;

    Info<< "\nDistributing faces to processors" << endl;

    // Loop through all internal faces and decide which processor they belong to
    // First visit all internal faces. If cells at both sides belong to the
    // same processor, the face is an internal face. If they are different,
    // it belongs to both processors.

    // procFaceAddressing_.setSize(nProcs_);

    // Internal faces
    forAll(neighbour, facei)
    {
        if (procNo_ == cellToProc_[owner[facei]])
        {
            if (cellToProc_[owner[facei]] == cellToProc_[neighbour[facei]])
            {
                // Face internal to processor. Notice no turning index.
                procFaceAddressing_.append(facei+1);
            }
        }
    }
    // for all processors, set the size of start index and patch size
    // lists to the number of patches in the mesh

    procPatchSize_.setSize(patches.size());
    procPatchStartIndex_.setSize(patches.size());

    forAll(patches, patchi)
    {
        // Reset size and start index for all processors

        procPatchSize_[patchi] = 0;
        procPatchStartIndex_[patchi] = procFaceAddressing_.size();

        const label patchStart = patches[patchi].start();

        if (!isA<cyclicPolyPatch>(patches[patchi]))
        {
            // Normal patch. Add faces to processor where the cell
            // next to the face lives

            const labelUList& patchFaceCells =
                patches[patchi].faceCells();

            forAll(patchFaceCells, facei)
            {
                const label curProc = cellToProc_[patchFaceCells[facei]];
                if (procNo_ == curProc)
                {
                    // add the face without turning index
                    procFaceAddressing_.append(patchStart+facei+1);

                    // increment the number of faces for this patch
                    procPatchSize_[patchi]++;
                }
            }
        }
        else
        {
            const cyclicPolyPatch& pp = refCast<const cyclicPolyPatch>
            (
                patches[patchi]
            );
            // cyclic: check opposite side on this processor
            const labelUList& patchFaceCells = pp.faceCells();

            const labelUList& nbrPatchFaceCells =
                pp.neighbPatch().faceCells();

            forAll(patchFaceCells, facei)
            {
                const label curProc = cellToProc_[patchFaceCells[facei]];
                const label nbrProc = cellToProc_[nbrPatchFaceCells[facei]];
                if (procNo_ == curProc)
                {
                    if (curProc == nbrProc)
                    {
                        // add the face without turning index
                        procFaceAddressing_.append(patchStart+facei+1);
                        // increment the number of faces for this patch
                        procPatchSize_[patchi]++;
                    }
                }
            }
        }
    }

    outputFilePtr.reset(new OFstream(path/"procPatchSize_"));
    outputFilePtr() << procPatchSize_ << endl;
    
    outputFilePtr.reset(new OFstream(path/"procPatchStartIndex_"));
    outputFilePtr() << procPatchStartIndex_ << endl;

    outputFilePtr.reset(new OFstream(path/"procFaceAddressing_"));
    outputFilePtr() << procFaceAddressing_ << endl;
    
}

bool Foam::parallelDomainDecomposition::writeDecomposition()
{
    Pout << "writeDecomposition is called" << endl;
    Pout << "\nCalculating original mesh data for process " << procNo_ << endl;

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

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}

/*
    TODO: there is a subtle bug that I could not fix yet.
    To understand the bug, run decomposePar and parallelDecompose.
    Then compare procCellAddressing_
*/
// assign each cell to its corresponding process.
// note: the code is copied from invertOneToMany method and changed.
Foam::labelList Foam::parallelDomainDecomposition::assignCellsToProc(
    const labelUList& map
)
{
    // each process will hold the information of how many cells that they have
    label procCellSize(0);

    for (const label newIdx : map)
    {
        if (newIdx == procNo_)
        {
            ++procCellSize;
        }
    }

    labelList inverse(procCellSize);
    procCellSize = 0;

    label i = 0;
    for (const label newIdx : map)
    {
        if (newIdx == procNo_)
        {
            inverse[procCellSize++] = i;
        }
        ++i;
    }

    return inverse;
}

// ************************************************************************* //
