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
        : public fvMesh
    {
    private:
        label procNo_;
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

        //- Optional: points at the facesInstance
        autoPtr<pointIOField> facesInstancePointsPtr_;
        DynamicList<label> procFaceAddressing_;

        void distributeCells();
        labelList assignCellsToProc(const labelUList& map);

    public:
        parallelDomainDecomposition(const IOobject& io, label procNo);
        ~parallelDomainDecomposition() = default;

        void decomposeMesh();
        bool writeDecomposition();

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
