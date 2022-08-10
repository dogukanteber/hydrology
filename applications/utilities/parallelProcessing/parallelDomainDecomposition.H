/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Class
    Foam::parallelDomainDecomposition

Description
    Automatic domain decomposition class for individual processes
    and finite-volume meshes

SourceFiles
    parallelDomainDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef parallelDomainDecomposition_H
#define parallelDomainDecomposition_H

#include "fvMesh.H"
#include "simpleGeomDecomp.H"
#include "labelList.H"
#include "point.H"
#include "Time.H"
#include "volFields.H"

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                         Class parallelDomainDecomposition Declaration
    \*---------------------------------------------------------------------------*/

    class parallelDomainDecomposition
        : public fvMesh,
          public simpleGeomDecomp
    {
    private:
        autoPtr<pointIOField> facesInstancePointsPtr_;
        label procNo_;
        label nProcs_;
        labelList cellToProc_;
        labelList procPointAddressing_;
        labelList procCellAddressing_;
        labelList procPatchSize_;
        labelList procPatchStartIndex_;
        labelList procNeighbourProcessors_;
        labelList procProcessorPatchSize_;
        labelList procProcessorPatchStartIndex_;
        labelListList procProcessorPatchSubPatchIDs_;
        labelListList procProcessorPatchSubPatchStarts_;

        DynamicList<label> procFaceAddressing_;

        void distributeCells();
        labelList assignCellsToProc(const labelUList& map);
        void addInterProcFace
        (
            const label facei,
            const label ownerProc,
            const label nbrProc,

            List<Map<label>>& nbrToInterPatch,
            List<DynamicList<DynamicList<label>>>& interPatchFaces
        ) const;

        //- Generate sub patch info for processor cyclics
        template<class BinaryOp>
        void processInterCyclics
        (
            const polyBoundaryMesh& patches,
            List<DynamicList<DynamicList<label>>>& interPatchFaces,
            List<Map<label>>& procNbrToInterPatch,
            labelListList& subPatchIDs,
            labelListList& subPatchStarts,
            bool owner,
            BinaryOp bop
        );

        void append(labelList& lst, const label elem);
        void mark
        (
            const labelList& zoneElems,
            const label zoneI,
            labelList& elementToZone
        );

        labelList simpleDecomposition(const pointField& cellPoints);
        std::pair<label, label> assignToProcessorGroup(labelList& processorGroup, const label nProcGroup);
        std::pair<label, label> nnn(labelList& processorGroup, const label nProcGroup);

    public:
        parallelDomainDecomposition(const IOobject& io, const IOdictionary& dict, label procNo, label nProcs);
        ~parallelDomainDecomposition() = default;

        void decomposeMesh();
        bool writeDecomposition();

        inline label nProcs() { return nProcs_; }
        inline label procNo() { return procNo_; }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
