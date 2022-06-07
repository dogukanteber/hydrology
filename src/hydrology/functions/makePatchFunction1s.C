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

#include "PatchFunction1.H"
#include "fieldTypes.H"
#include "UniformValueField.H"
#include "addToRunTimeSelectionTable.H"

#if (OPENFOAM < 2112)
// Migrated to UniformValueField.H for OpenFOAM-v2112
#define addUniformValueFieldFunction1s(F1Type, Type)                           \
    PatchFunction1<Type>::adddictionaryConstructorToTable                      \
    <PatchFunction1Types::UniformValueField<Type>>                             \
        add##F1Type##UniformValueField##Type##ConstructorToTable_(#F1Type);
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Add UniformFieldValue under the same name as Function1
    // See makeFunction1s.C

    addUniformValueFieldFunction1s(sine1, scalar);
    addUniformValueFieldFunction1s(sine1, vector);
    addUniformValueFieldFunction1s(sine1, sphericalTensor);
    addUniformValueFieldFunction1s(sine1, symmTensor);
    addUniformValueFieldFunction1s(sine1, tensor);

    addUniformValueFieldFunction1s(cosine1, scalar);
    addUniformValueFieldFunction1s(cosine1, vector);
    addUniformValueFieldFunction1s(cosine1, sphericalTensor);
    addUniformValueFieldFunction1s(cosine1, symmTensor);
    addUniformValueFieldFunction1s(cosine1, tensor);

    addUniformValueFieldFunction1s(square1, scalar);
    addUniformValueFieldFunction1s(square1, vector);
    addUniformValueFieldFunction1s(square1, sphericalTensor);
    addUniformValueFieldFunction1s(square1, symmTensor);
    addUniformValueFieldFunction1s(square1, tensor);

    addUniformValueFieldFunction1s(seasonal, scalar);
    addUniformValueFieldFunction1s(seasonal, vector);
    addUniformValueFieldFunction1s(seasonal, sphericalTensor);
    addUniformValueFieldFunction1s(seasonal, symmTensor);
    addUniformValueFieldFunction1s(seasonal, tensor);

}


// ************************************************************************* //
