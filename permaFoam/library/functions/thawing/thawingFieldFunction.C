/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "thawingFieldFunction.H"
#include "SubField.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"


// As per OpenFOAM-v2012

#ifndef makeConcretePatchFunction1Type
// Define (non-templated) PatchFunction1, add to (templated) run-time selection
#define makeConcretePatchFunction1Type(SS, Type)                               \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    PatchFunction1<Type>::adddictionaryConstructorToTable                      \
        <PatchFunction1Types::SS>                                              \
        add##SS##Type##ConstructorToTable_;
#endif

#ifndef makeScalarPatchFunction1
// Define scalar PatchFunction1 and add to (templated) run-time selection
#define makeScalarPatchFunction1(SS)                                           \
    makeConcretePatchFunction1Type(SS, scalar);
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

makeScalarPatchFunction1(thawingFieldFunction);

} // End namespace PatchFunction1Types
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PatchFunction1Types::thawingFieldFunction::thawingFieldFunction
(
    const polyPatch& pp,
    const word& type,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<scalar>(pp, entryName, dict, faceValues),
    faceName_(dict.get<word>("Tface")),
    meltName_(),
    meltValue_(0),
    shiftValue_(dict.getOrDefault<scalar>("Tshift", 0))
{
    token tok(dict.lookup("Tmelt"));

    if (tok.isNumber())
    {
        dict.readEntry("Tmelt", meltValue_);
    }
    else
    {
        dict.readEntry("Tmelt", meltName_);
    }
}


Foam::PatchFunction1Types::thawingFieldFunction::thawingFieldFunction
(
    const thawingFieldFunction& rhs
)
:
    PatchFunction1<scalar>(rhs),
    faceName_(rhs.faceName_),
    meltName_(rhs.meltName_),
    meltValue_(rhs.meltValue_),
    shiftValue_(rhs.shiftValue_)
{}


Foam::PatchFunction1Types::thawingFieldFunction::thawingFieldFunction
(
    const thawingFieldFunction& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<scalar>(rhs),
    faceName_(rhs.faceName_),
    meltName_(rhs.meltName_),
    meltValue_(rhs.meltValue_),
    shiftValue_(rhs.shiftValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::PatchFunction1Types::thawingFieldFunction::value
(
    const scalar x
) const
{
    // ONLY valid for faceValues !!!

    const label patchi = this->patch_.index();

    const fvMesh* fvmesh = isA<fvMesh>(this->patch_.boundaryMesh().mesh());

    if (!fvmesh)
    {
        FatalErrorInFunction
            << "Not an fvMesh!" << nl
            << exit(FatalError);
    }


    // Tface
    const auto* sfield = fvmesh->cfindObject<surfaceScalarField>(faceName_);

    if (!sfield)
    {
        FatalErrorInFunction
            << "No Tface field: " << faceName_ << nl
            << exit(FatalError);
    }

    const fvsPatchField<scalar>& T_faces = sfield->boundaryField()[patchi];


    // Tmelt
    const auto* vfield = fvmesh->cfindObject<volScalarField>(meltName_);

    if (vfield)
    {
        const fvPatchField<scalar>& T_melt = vfield->boundaryField()[patchi];

        return pos0(T_faces - (T_melt + shiftValue_));
    }
    else
    {
        if (!meltName_.empty())
        {
            FatalErrorInFunction
                << "No Tmelt field: " << meltName_ << nl
                << exit(FatalError);
        }
    }


    return pos0(T_faces - (meltValue_ + shiftValue_));
}


void Foam::PatchFunction1Types::thawingFieldFunction::autoMap
(
    const FieldMapper& mapper
)
{
    PatchFunction1<scalar>::autoMap(mapper);
}


void Foam::PatchFunction1Types::thawingFieldFunction::rmap
(
    const PatchFunction1<scalar>& pf1,
    const labelList& addr
)
{
    PatchFunction1<scalar>::rmap(pf1, addr);
}


void Foam::PatchFunction1Types::thawingFieldFunction::writeData
(
    Ostream& os
) const
{
    os.writeEntry("Tface", faceName_);

    if (meltName_.empty())
    {
        // Value
        os.writeEntry("Tmelt", meltValue_);
    }
    else
    {
        // Field name
        os.writeEntry("Tmelt", meltName_);
    }

    os.writeEntryIfDifferent<scalar>("Tshift", shiftValue_, 0);
}


// ************************************************************************* //
