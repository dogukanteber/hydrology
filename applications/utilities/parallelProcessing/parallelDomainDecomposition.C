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
#include "cyclicPolyPatch.H"

// change the constructor
Foam::parallelDomainDecomposition::parallelDomainDecomposition(
    const IOobject &io,
    label procNo,
    label nProcs)
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
      nProcs_(nProcs)
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
        label owner = cellToProc_[owner[facei]];
        if (procNo_ == owner)
        {
            if (owner == cellToProc_[neighbour[facei]])
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

    // Done internal bits of the new mesh and the ordinary patches.

    // Below two lists are serial!!
    // Per processor, from neighbour processor to the inter-processor patch
    // that communicates with that neighbour
    List<Map<label>> procNbrToInterPatch(nProcs_);

    // Per processor the faces per inter-processor patch
    List<DynamicList<DynamicList<label>>> interPatchFaces(nProcs_);

    // Processor boundaries from internal faces
    forAll(neighbour, facei)
    {
        label ownerProc = cellToProc_[owner[facei]];
        label nbrProc = cellToProc_[neighbour[facei]];

        if (ownerProc == procNo_)
        {
            if (ownerProc != nbrProc)
            {
                // inter - processor patch face found.
                addInterProcFace
                (
                    facei,
                    ownerProc,
                    nbrProc,

                    procNbrToInterPatch,
                    interPatchFaces
                );
            }
        }
    }

    outputFilePtr.reset(new OFstream(outputDir/"procNbrToInterPatch"));
    outputFilePtr() << procNbrToInterPatch << endl;

    outputFilePtr.reset(new OFstream(outputDir/"interPatchFaces"));
    outputFilePtr() << interPatchFaces << endl;

    // Add the proper processor faces to the sub information. For faces
    // originating from internal faces this is always -1.
    labelListList subPatchIDs;
    labelListList subPatchStarts;

    label nInterfaces = interPatchFaces[procNo_].size();

    subPatchIDs.setSize(nInterfaces, labelList(1, label(-1)));
    subPatchStarts.setSize(nInterfaces, labelList(1, Zero));

    outputFilePtr.reset(new OFstream(outputDir/"subPatchIDs"));
    outputFilePtr() << subPatchIDs << endl;

    outputFilePtr.reset(new OFstream(outputDir/"subPatchStarts"));
    outputFilePtr() << subPatchStarts << endl;

    // Special handling needed for the case that multiple processor cyclic
    // patches are created on each local processor domain, e.g. if a 3x3 case
    // is decomposed using the decomposition:
    //
    //              | 1 | 0 | 2 |
    //  cyclic left | 2 | 0 | 1 | cyclic right
    //              | 2 | 0 | 1 |
    //
    // - processors 1 and 2 will both have pieces of both cyclic left- and
    //   right sub-patches present
    // - the interface patch faces are stored in a single list, where each
    //   sub-patch is referenced into the list using a patch start index and
    //   size
    // - if the patches are in order (in the boundary file) of left, right
    //   - processor 1 will send: left, right
    //   - processor 1 will need to receive in reverse order: right, left
    //   - similarly for processor 2
    // - the sub-patches are therefore generated in 4 passes of the patch lists
    //   1. add faces from owner patch where local proc i < nbr proc i
    //   2. add faces from nbr patch where local proc i < nbr proc i
    //   3. add faces from owner patch where local proc i > nbr proc i
    //   4. add faces from nbr patch where local proc i > nbr proc i

    processInterCyclics
    (
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        true,
        lessOp<label>()
    );

    processInterCyclics
    (
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        false,
        lessOp<label>()
    );

    processInterCyclics
    (
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        false,
        greaterOp<label>()
    );

    processInterCyclics
    (
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        true,
        greaterOp<label>()
    );

    outputFilePtr.reset(new OFstream(outputDir/"procNbrToInterPatch"));
    outputFilePtr() << procNbrToInterPatch << endl;
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
    label len(4);
    labelList sizes(len, Zero);

    for (const label newIdx : map)
    {
        if (newIdx >= 0)
        {
            ++sizes[newIdx];
        }
    }

    Info << sizes[0] << tab << sizes[1] << tab << sizes[2] << tab << sizes[3] << endl;
    Pout << "P: " << sizes[0] << tab << sizes[1] << tab << sizes[2] << tab << sizes[3] << endl;
    // autoPtr<OFstream> outputFilePtr;
    // fileName path("/home/dogukan/Documents/hydrology/tutorials/permaFoam/demoCase/debugprocessor" + Foam::name(procNo_));
    // outputFilePtr.reset(new OFstream(path/"assingCellsToProc"));
    // each process will hold the information of how many cells that they have
    label procCellSize(0);
    // label count(1);
    for (const label newIdx : map)
    {
        // outputFilePtr() << "count : " << count << tab << "newIdx : " << newIdx << tab;
        if (newIdx == procNo_)
        {
            ++procCellSize;
            // outputFilePtr() << "Inside the if. newIdx : " << newIdx << " procNo : " << procNo_ << " procCellSize : " << procCellSize << tab << endl << endl;
        }
        // count++;
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

void Foam::parallelDomainDecomposition::addInterProcFace
(
    const label facei,
    const label ownerProc,
    const label nbrProc,

    List<Map<label>>& nbrToInterPatch,
    List<DynamicList<DynamicList<label>>>& interPatchFaces
) const
{
    // Introduce turning index only for internal faces (are duplicated).
    const label ownerIndex = facei+1;
    const label nbrIndex = -(facei+1);

    const auto patchiter = nbrToInterPatch[ownerProc].cfind(nbrProc);

    if (patchiter.found())
    {
        // Existing interproc patch. Add to both sides.
        const label toNbrProcPatchi = *patchiter;
        interPatchFaces[ownerProc][toNbrProcPatchi].append(ownerIndex);

        if (isInternalFace(facei))
        {
            label toOwnerProcPatchi = nbrToInterPatch[nbrProc][ownerProc];
            interPatchFaces[nbrProc][toOwnerProcPatchi].append(nbrIndex);
        }
    }
    else
    {
        // Create new interproc patches.
        const label toNbrProcPatchi = nbrToInterPatch[ownerProc].size();
        nbrToInterPatch[ownerProc].insert(nbrProc, toNbrProcPatchi);

        DynamicList<label> oneFace;
        oneFace.append(ownerIndex);
        interPatchFaces[ownerProc].append(oneFace);

        if (isInternalFace(facei))
        {
            label toOwnerProcPatchi = nbrToInterPatch[nbrProc].size();
            nbrToInterPatch[nbrProc].insert(ownerProc, toOwnerProcPatchi);
            oneFace.clear();
            oneFace.append(nbrIndex);
            interPatchFaces[nbrProc].append(oneFace);
        }
    }
}

template<class BinaryOp>
void Foam::parallelDomainDecomposition::processInterCyclics
(
    const polyBoundaryMesh& patches,
    List<DynamicList<DynamicList<label>>>& interPatchFaces,
    List<Map<label>>& procNbrToInterPatch,
    labelListList& subPatchIDs,
    labelListList& subPatchStarts,
    bool owner,
    BinaryOp bop
) const
{
    // Processor boundaries from split cyclics
    forAll(patches, patchi)
    {
        const auto& pp = patches[patchi];
        const auto* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner() == owner)
        {
            // cyclic: check opposite side on this processor
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            // cyclic: check opposite side on this processor
            const labelUList& patchFaceCells = pp.faceCells();
            const labelUList& nbrPatchFaceCells = nbrPatch.faceCells();

            // Store old sizes. Used to detect which inter-proc patches
            // have been added to.
            labelList oldInterfaceSizes;
            labelList& curOldSizes = oldInterfaceSizes;

            curOldSizes.setSize(interPatchFaces[procNo_].size());
            forAll(curOldSizes, interI)
            {
                curOldSizes[interI] =
                    interPatchFaces[procNo_][interI].size();
            }

            // Add faces with different owner and neighbour processors
            forAll(patchFaceCells, facei)
            {
                const label ownerProc = cellToProc_[patchFaceCells[facei]];
                const label nbrProc = cellToProc_[nbrPatchFaceCells[facei]];
                if (bop(ownerProc, nbrProc))
                {
                    // inter - processor patch face found.
                    addInterProcFace
                    (
                        pp.start()+facei,
                        ownerProc,
                        nbrProc,
                        procNbrToInterPatch,
                        interPatchFaces
                    );
                }
            }

            // 1. Check if any faces added to existing interfaces
            const labelList& curOldSizes = oldInterfaceSizes;

            forAll(curOldSizes, interI)
            {
                label oldSz = curOldSizes[interI];
                if (interPatchFaces[procNo_][interI].size() > oldSz)
                {
                    // Added faces to this interface. Add an entry
                    append(subPatchIDs[interI], patchi);
                    append(subPatchStarts[interI], oldSz);
                }
            }
            // 2. Any new interfaces
            label nIntfcs = interPatchFaces[procNo_].size();
            subPatchIDs.setSize(nIntfcs, labelList(1, patchi));
            subPatchStarts.setSize(nIntfcs, labelList(1, Zero));
        }
    }
}

// ************************************************************************* //
