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
#include "IOstreams.H"
#include "bitSet.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "OSspecific.H"
#include "Map.H"
#include "DynamicList.H"
#include "fvFieldDecomposer.H"
#include "IOobjectList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "OPstream.H"
#include "IPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parallelDomainDecomposition::parallelDomainDecomposition(
    const IOobject &io,
    const IOdictionary& dict,
    label procNo,
    label nProcs)
    : fvMesh(io),
      simpleGeomDecomp(dict),
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
      nProcs_(nProcs),
      cellToProc_(4250)
{
    // Pout << "parallelDomainDecomposition constructor is called for process " << procNo_ << endl;
}

// A list compare binary predicate for normal sort by vector component
struct vectorLessOp
{
    const UList<vector>& values;
    direction sortCmpt;

    vectorLessOp(const UList<vector>& list, direction cmpt = vector::X)
    :
        values(list),
        sortCmpt(cmpt)
    {}

    void setComponent(direction cmpt)
    {
        sortCmpt = cmpt;
    }

    bool operator()(const label a, const label b) const
    {
        return values[a][sortCmpt] < values[b][sortCmpt];
    }
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*
    !!! CAUTION !!!
    This method is tested only when nProcGroup == nProcs and nProcGroup == 1
    It is left this way because the time strict
*/
std::pair<label, label> Foam::parallelDomainDecomposition::assignToProcessorGroup(
    labelList& processorGroup,
    const label nProcGroup
)
{
    const label nProcs = Pstream::nProcs();
    const label procNo = Pstream::myProcNo();
    const label totalCells = processorGroup.size();

    std::pair<label, label> range;
    label startIndex = 0, endIndex = 0;
    if (nProcs == nProcGroup)
    {
        // distribute cells evenly among each process
        label cellsPerProc = totalCells / nProcs;
        label remainder = totalCells % nProcs;
        if (procNo < remainder)
        {
            startIndex = procNo * (cellsPerProc + 1);
            endIndex = startIndex + (cellsPerProc + 1);
        }
        else
        {
            startIndex = procNo * cellsPerProc + remainder;
            endIndex = startIndex + cellsPerProc;
        }
        range.first = startIndex;
        range.second = endIndex;
        // Pout << startIndex << " - " << endIndex << endl;
        for (; startIndex < endIndex; startIndex++)
        {
            processorGroup[startIndex] = procNo;
        }
    }
    else
    {
        // when nProcGroups is 1, distribute the workload among processes
        label cellsPerProc = totalCells / nProcs;
        label remainder = totalCells % nProcs;
        if (procNo < remainder)
        {
            startIndex = procNo * (cellsPerProc + 1);
            endIndex = startIndex + (cellsPerProc + 1);
        }
        else
        {
            startIndex = procNo * cellsPerProc + remainder;
            endIndex = startIndex + cellsPerProc;
        }
        range.first = startIndex;
        range.second = endIndex;
        for (; startIndex < endIndex; startIndex++)
        {
            processorGroup[startIndex] = 0;
        }
    }

    return range;
}

Foam::labelList Foam::parallelDomainDecomposition::simpleDecomposition(const pointField& cellPoints)
{
    // fileName outputDir = cwd()/"simpleDecomp";
    // mkDir(outputDir);
    // autoPtr<OFstream> outputFilePtr;

    labelList finalDecomp(cellPoints.size(), -1);

    labelList processorGroups(cellPoints.size());

    labelList pointIndices(identity(cellPoints.size()));

    const pointField rotatedPoints(adjustPoints(cellPoints));

    vectorLessOp sorter(rotatedPoints);

    sorter.setComponent(vector::X);
    Foam::sort(pointIndices, sorter);

    std::pair<label, label> range = assignToProcessorGroup(processorGroups, n_.x());

    for (label i=range.first; i<range.second; i++)
    {
        finalDecomp[pointIndices[i]] = processorGroups[i];
    }

    sorter.setComponent(vector::Y);
    Foam::sort(pointIndices, sorter);

    range = assignToProcessorGroup(processorGroups, n_.y());

    for (label i=range.first; i<range.second; i++)
    {
        finalDecomp[pointIndices[i]] += n_.x()*processorGroups[i];
    }

    sorter.setComponent(vector::Z);
    Foam::sort(pointIndices, sorter);

    range = assignToProcessorGroup(processorGroups, n_.z());

    for (label i=range.first; i<range.second; i++)
    {
        finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
    }

    labelListList n_finalDecomp(UPstream::nProcs());
    n_finalDecomp[UPstream::myProcNo()] = finalDecomp;
    Pstream::gatherList(n_finalDecomp);

    labelList combined(
        ListListOps::combine<labelList>(n_finalDecomp, accessOp<labelList>())
    );

    labelList finalVersion(cellPoints.size());
    if (Pstream::master())
    {
        const label meshSize = cellPoints.size();
        const label totalMeshSize = cellPoints.size() * Pstream::nProcs();

        label start = 0;
        label end = totalMeshSize;
        label index = 0;
        while ( start < end )
        {
            label currentElement = combined[start];
            label invalid(-1);
            while (currentElement != invalid)
            {
                finalVersion[index] = currentElement;
                index++;
                start++;
                // Info << start << endl;
                if (start == end)
                {
                    break;
                }
                currentElement = combined[start];
            }
            if (start == end)
            {
                break;
            }
            if (finalVersion[index-1] == Pstream::lastSlave())
            {
                start = start - (meshSize * (Pstream::nProcs() - 1));
            }
            else
            {
                start += meshSize;
            }
        }
    }
    Pstream::scatter(finalVersion);
    
    // outputFilePtr.reset(new OFstream(outputDir/"finalVersion"));
    // outputFilePtr() << finalVersion << endl;

    return finalVersion;
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

    cellToProc_ = simpleDecomposition(this->cellCentres());

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
    // autoPtr<OFstream> outputFilePtr;
    // fileName path("/home/dogukan/Documents/hydrology/tutorials/permaFoam/demoCase/debugprocessor" + Foam::name(procNo_));
    // mkDir(path);
    // outputFilePtr.reset(new OFstream(path/"cellToProc_"));
    // outputFilePtr() << map << endl;

    // cellToProc_ comes the same for each process. Tested and verified.

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
)
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
                if (ownerProc == procNo_)
                {
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
            }

            // 1. Check if any faces added to existing interfaces
            const labelList& curOldSizes2 = oldInterfaceSizes;

            forAll(curOldSizes2, interI)
            {
                label oldSz = curOldSizes2[interI];
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

void Foam::parallelDomainDecomposition::append(labelList& lst, const label elem)
{
    label sz = lst.size();
    lst.setSize(sz+1);
    lst[sz] = elem;
}

void Foam::parallelDomainDecomposition::mark
(
    const labelList& zoneElems,
    const label zoneI,
    labelList& elementToZone
)
{
    for (const label pointi : zoneElems)
    {
        if (elementToZone[pointi] == -1)
        {
            // First occurrence
            elementToZone[pointi] = zoneI;
        }
        else if (elementToZone[pointi] >= 0)
        {
            // Multiple zones
            elementToZone[pointi] = -2;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parallelDomainDecomposition::decomposeMesh()
{
    fileName path = cwd()/"simpleDecomp";
    mkDir(path);
    autoPtr<OFstream> outputFilePtr;

    distributeCells();

    outputFilePtr.reset(new OFstream(path/"cellToProc_"));
    outputFilePtr() << cellToProc_ << endl;

    // cellToProc_ is tested and verified. Yay!

    Info<< "\nCalculating original mesh data" << endl;

    // set references to the original mesh
    const polyBoundaryMesh& patches = boundaryMesh();
    const faceList& fcs = faces();
    const labelList& owner = faceOwner();
    const labelList& neighbour = faceNeighbour();

    Info<< "\nDistributing cells to processors" << endl;
    procCellAddressing_ = assignCellsToProc(cellToProc_);

    // procCellAddressing_ is working. Tested and verified.

    // test code
    // labelListList combined(Pstream::nProcs());
    // combined[Pstream::myProcNo()] = procCellAddressing_;

    // Pstream::gatherList(combined);

    // if (Pstream::master())
    // {
    //     outputFilePtr.reset(new OFstream(path/"procCellAddressing_"));
    //     outputFilePtr() << combined << endl;
    // }

    Info<< "\nDistributing faces to processors" << endl;

    // Loop through all internal faces and decide which processor they belong to
    // First visit all internal faces. If cells at both sides belong to the
    // same processor, the face is an internal face. If they are different,
    // it belongs to both processors.

    // procFaceAddressing_.setSize(nProcs_);

    // Internal faces
    forAll(neighbour, facei)
    {
        label currProc = cellToProc_[owner[facei]];
        if (procNo_ == currProc)
        {
            if (currProc == cellToProc_[neighbour[facei]])
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

    // BELOW ARE TESTED AND VERIFIED
    /*
        procFaceAddressing_
        procPatchStartIndex_
        procFaceAddressing_
    */

    // Done internal bits of the new mesh and the ordinary patches.

    // Below two lists are serial!!
    // Per processor, from neighbour processor to the inter-processor patch
    // that communicates with that neighbour
    List<Map<label>> procNbrToInterPatch(nProcs_);

    // Per processor the faces per inter-processor patch
    List<DynamicList<DynamicList<label>>> interPatchFaces(nProcs_);

    // below can be improved in terms of perfomance
    if (Pstream::master())
    {
        // Processor boundaries from internal faces
        forAll(neighbour, facei)
        {
            label ownerProc = cellToProc_[owner[facei]];
            label nbrProc = cellToProc_[neighbour[facei]];

            // if (ownerProc == procNo_)
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
    }
    Pstream::scatter(procNbrToInterPatch);
    Pstream::scatter(interPatchFaces);

    // procNbrToInterPatch and interPatchFaces is working correctly!

    // Add the proper processor faces to the sub information. For faces
    // originating from internal faces this is always -1.
    labelListList subPatchIDs;
    labelListList subPatchStarts;

    label nInterfaces = interPatchFaces[procNo_].size();

    subPatchIDs.setSize(nInterfaces, labelList(1, label(-1)));
    subPatchStarts.setSize(nInterfaces, labelList(1, Zero));

    // subPatchIDs and subPatchStarts are working correctly!

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

    // I think procNbrToInterPatch is working correctly. So, I did not test it.

    // Sort inter-proc patch by neighbour
    labelList order;
    label nInterfaces2 = procNbrToInterPatch[procNo_].size();

    procNeighbourProcessors_.setSize(nInterfaces2);
    procProcessorPatchSize_.setSize(nInterfaces2);
    procProcessorPatchStartIndex_.setSize(nInterfaces2);
    procProcessorPatchSubPatchIDs_.setSize(nInterfaces2);
    procProcessorPatchSubPatchStarts_.setSize(nInterfaces2);

    // Get sorted neighbour processors
    const Map<label>& curNbrToInterPatch = procNbrToInterPatch[procNo_];
    labelList nbrs = curNbrToInterPatch.toc();

    sortedOrder(nbrs, order);

    DynamicList<DynamicList<label>>& curInterPatchFaces =
        interPatchFaces[procNo_];

    forAll(nbrs, i)
    {
        const label nbrProc = nbrs[i];
        const label interPatch = curNbrToInterPatch[nbrProc];

        procNeighbourProcessors_[i] = nbrProc;
        procProcessorPatchSize_[i] =
            curInterPatchFaces[interPatch].size();
        procProcessorPatchStartIndex_[i] =
            procFaceAddressing_.size();

        // Add size as last element to substarts and transfer
        append
        (
            subPatchStarts[interPatch],
            curInterPatchFaces[interPatch].size()
        );
        procProcessorPatchSubPatchIDs_[i].transfer
        (
            subPatchIDs[interPatch]
        );
        procProcessorPatchSubPatchStarts_[i].transfer
        (
            subPatchStarts[interPatch]
        );

        DynamicList<label>& interPatchFaces =
            curInterPatchFaces[interPatch];

        forAll(interPatchFaces, j)
        {
            procFaceAddressing_.append(interPatchFaces[j]);
        }
        interPatchFaces.clearStorage();
    }
    curInterPatchFaces.clearStorage();
    procFaceAddressing_.shrink();

    // below are working correctly:

    /*
        procNeighbourProcessors_
        procProcessorPatchSize_
        procNbrToInterPatch
        procProcessorPatchSubPatchIDs_
        procProcessorPatchSubPatchStarts_  
    */

    Info<< "\nDistributing points to processors" << endl;
    // For every processor, loop through the list of faces for the processor.
    // For every face, loop through the list of points and mark the point as
    // used for the processor. Collect the list of used points for the
    // processor.

    bitSet pointsInUse(nPoints(), false);

    // For each of the faces used
    for (const label facei : procFaceAddressing_)
    {
        // Because of the turning index, some face labels may be -ve
        const labelList& facePoints = fcs[mag(facei) - 1];

        // Mark the face points as being used
        pointsInUse.set(facePoints);
    }
    procPointAddressing_ = pointsInUse.sortedToc();

    // procPointAddressing_ is working as well.
}

bool Foam::parallelDomainDecomposition::writeDecomposition()
{
    Info<< "\nConstructing processor meshes" << endl;

    // Mark point/faces/cells that are in zones.
    // -1   : not in zone
    // -2   : in multiple zones
    // >= 0 : in single given zone
    // This will give direct lookup of elements that are in a single zone
    // and we'll only have to revert back to searching through all zones
    // for the duplicate elements

    // contains -1
    // Point zones
    labelList pointToZone(points().size(), -1);

    forAll(pointZones(), zonei)
    {
        mark(pointZones()[zonei], zonei, pointToZone);
    }

    // contains -1
    // Face zones
    labelList faceToZone(faces().size(), -1);

    forAll(faceZones(), zonei)
    {
        mark(faceZones()[zonei], zonei, faceToZone);
    }

    // contains -1
    // Cell zones
    labelList cellToZone(nCells(), -1);

    forAll(cellZones(), zonei)
    {
        mark(cellZones()[zonei], zonei, cellToZone);
    }

    PtrList<const cellSet> cellSets;
    PtrList<const faceSet> faceSets;
    PtrList<const pointSet> pointSets;

    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    // Write out the meshes
    
    // holds reference to processorN/constant/polyMesh/points
    // Create processor points
    const labelList& curPointLabels = procPointAddressing_;

    const pointField& meshPoints = points();

    labelList pointLookup(nPoints(), -1);

    pointField procPoints(curPointLabels.size());

    forAll(curPointLabels, pointi)
    {
        procPoints[pointi] = meshPoints[curPointLabels[pointi]];

        pointLookup[curPointLabels[pointi]] = pointi;
    }

    // Create processor faces
    const labelList& curFaceLabels = procFaceAddressing_;

    const faceList& meshFaces = faces();

    // labelList faceLookup(nFaces(), -1);
    labelList faceLookup(nFaces(), -1);

    faceList procFaces(curFaceLabels.size());

    forAll(curFaceLabels, facei)
    {
        // Mark the original face as used
        // Remember to decrement the index by one (turning index)
        label curF = mag(curFaceLabels[facei]) - 1;

        faceLookup[curF] = facei;

        // get the original face
        labelList origFaceLabels;

        if (curFaceLabels[facei] >= 0)
        {
            // face not turned
            origFaceLabels = meshFaces[curF];
        }
        else
        {
            origFaceLabels = meshFaces[curF].reverseFace();
        }

        // translate face labels into local point list
        face& procFaceLabels = procFaces[facei];

        procFaceLabels.setSize(origFaceLabels.size());

        forAll(origFaceLabels, pointi)
        {
            procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
        }
    }

    // Create processor cells
    const labelList& curCellLabels = procCellAddressing_;

    const cellList& meshCells = cells();

    cellList procCells(curCellLabels.size());

    forAll(curCellLabels, celli)
    {
        const labelList& origCellLabels = meshCells[curCellLabels[celli]];

        cell& curCell = procCells[celli];

        curCell.setSize(origCellLabels.size());

        forAll(origCellLabels, cellFacei)
        {
            curCell[cellFacei] = faceLookup[origCellLabels[cellFacei]];
        }
    }

    // Create processor mesh without a boundary
    fileName processorCasePath
    (
        time().caseName()/("processor" + Foam::name(procNo_))
    );

    // create a database
    Time processorDb
    (
        Time::controlDictName,
        time().rootPath(),
        processorCasePath,
        word("system"),
        word("constant"),
        false
    );
    processorDb.setTime(time());

    // create the mesh. Two situations:
    // - points and faces come from the same time ('instance'). The mesh
    //   will get constructed in the same instance.
    // - points come from a different time (moving mesh cases).
    //   It will read the points belonging to the faces instance and
    //   construct the procMesh with it which then gets handled as above.
    //   (so with 'old' geometry).
    //   Only at writing time will it additionally write the current
    //   points.

    autoPtr<polyMesh> procMeshPtr;

    if (facesInstancePointsPtr_)
    {
        // Construct mesh from facesInstance.
        pointField facesInstancePoints
        (
            facesInstancePointsPtr_(),
            curPointLabels
        );

        procMeshPtr = autoPtr<polyMesh>::New
        (
            IOobject
            (
                this->polyMesh::name(), // region of undecomposed mesh
                facesInstance(),
                processorDb,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            std::move(facesInstancePoints),
            std::move(procFaces),
            std::move(procCells)
        );
    }
    else
    {
        // Pout << procPoints << endl;
        // Pout << procFaces << endl;
        // Pout << procCells << endl;
        procMeshPtr = autoPtr<polyMesh>::New
        (
            IOobject
            (
                this->polyMesh::name(), // region of undecomposed mesh
                facesInstance(),
                processorDb,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            std::move(procPoints),
            std::move(procFaces),
            // procCells variable contains negative values and
            // because of that the program is throwing an error
            std::move(procCells)
        );
    }
    polyMesh& procMesh = procMeshPtr();

    // Create processor boundary patches
    const labelList& curPatchSizes = procPatchSize_;

    const labelList& curPatchStarts = procPatchStartIndex_;

    const labelList& curNeighbourProcessors =
        procNeighbourProcessors_;

    const labelList& curProcessorPatchSizes =
        procProcessorPatchSize_;

    const labelList& curProcessorPatchStarts =
        procProcessorPatchStartIndex_;

    const labelListList& curSubPatchIDs =
        procProcessorPatchSubPatchIDs_;

    const labelListList& curSubStarts =
        procProcessorPatchSubPatchStarts_;

    const polyPatchList& meshPatches = boundaryMesh();

    // Count the number of inter-proc patches
    label nInterProcPatches = 0;
    forAll(curSubPatchIDs, procPatchi)
    {
        nInterProcPatches += curSubPatchIDs[procPatchi].size();
    }

    PtrList<polyPatch> procPatches
    (
        curPatchSizes.size() + nInterProcPatches
    );

    label nPatches = 0;

    forAll(curPatchSizes, patchi)
    {
        // Get the face labels consistent with the field mapping
        // (reuse the patch field mappers)
        const polyPatch& meshPatch = meshPatches[patchi];

        fvFieldDecomposer::patchFieldDecomposer patchMapper
        (
            SubList<label>
            (
                curFaceLabels,
                curPatchSizes[patchi],
                curPatchStarts[patchi]
            ),
            meshPatch.start()
        );

        // Map existing patches
        procPatches.set
        (
            nPatches,
            meshPatch.clone
            (
                procMesh.boundaryMesh(),
                nPatches,
                patchMapper.directAddressing(),
                curPatchStarts[patchi]
            )
        );

        nPatches++;
    }

    forAll(curProcessorPatchSizes, procPatchi)
    {
        const labelList& subPatchID = curSubPatchIDs[procPatchi];
        const labelList& subStarts = curSubStarts[procPatchi];

        label curStart = curProcessorPatchStarts[procPatchi];

        forAll(subPatchID, i)
        {
            label size =
            (
                i < subPatchID.size()-1
                ? subStarts[i+1] - subStarts[i]
                : curProcessorPatchSizes[procPatchi] - subStarts[i]
            );

            if (subPatchID[i] == -1)
            {
                // From internal faces
                procPatches.set
                (
                    nPatches,
                    new processorPolyPatch
                    (
                        size,
                        curStart,
                        nPatches,
                        procMesh.boundaryMesh(),
                        procNo_,
                        curNeighbourProcessors[procPatchi]
                    )
                );
            }
            else
            {
                const coupledPolyPatch& pcPatch
                    = refCast<const coupledPolyPatch>
                        (
                            boundaryMesh()[subPatchID[i]]
                        );

                procPatches.set
                (
                    nPatches,
                    new processorCyclicPolyPatch
                    (
                        size,
                        curStart,
                        nPatches,
                        procMesh.boundaryMesh(),
                        procNo_,
                        curNeighbourProcessors[procPatchi],
                        pcPatch.name(),
                        pcPatch.transform()
                    )
                );
            }

            curStart += size;
            ++nPatches;
        }
    }
    // Pout << procPatches << endl;
    // Info << procMesh << endl;
    procMesh.addPatches(procPatches);

    // Point zones
    {
        // const pointZoneMesh& pz = pointZones();

        // Go through all the zoned points and find out if they
        // belong to a zone.  If so, add it to the zone as
        // necessary
        // List<DynamicList<label>> zonePoints(pz.size());

        // // Estimate size
        // forAll(zonePoints, zonei)
        // {
        //     zonePoints[zonei].setCapacity(pz[zonei].size()/nProcs_);
        // }

    //     // Use the pointToZone map to find out the single zone (if any),
    //     // use slow search only for shared points.
    //     forAll(curPointLabels, pointi)
    //     {
    //         label curPoint = curPointLabels[pointi];

    //         label zonei = pointToZone[curPoint];

    //         if (zonei >= 0)
    //         {
    //             // Single zone.
    //             zonePoints[zonei].append(pointi);
    //         }
    //         else if (zonei == -2)
    //         {
    //             // Multiple zones. Lookup.
    //             forAll(pz, zonei)
    //             {
    //                 label index = pz[zonei].whichPoint(curPoint);

    //                 if (index != -1)
    //                 {
    //                     zonePoints[zonei].append(pointi);
    //                 }
    //             }
    //         }
    //     }

    //     procMesh.pointZones().clearAddressing();
    //     procMesh.pointZones().setSize(zonePoints.size());
    //     forAll(zonePoints, zonei)
    //     {
    //         procMesh.pointZones().set
    //         (
    //             zonei,
    //             pz[zonei].clone
    //             (
    //                 procMesh.pointZones(),
    //                 zonei,
    //                 zonePoints[zonei].shrink()
    //             )
    //         );
    //     }

    //     if (pz.size())
    //     {
    //         // Force writing on all processors
    //         procMesh.pointZones().writeOpt(IOobject::AUTO_WRITE);
    //     }
    }

    

    return true;
}


// ************************************************************************* //
