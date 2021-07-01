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

#include "Sine1.H"
#include "Cosine1.H"
#include "Square1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

makeFunction1Type(Cosine1, scalar);
makeFunction1Type(Cosine1, vector);
makeFunction1Type(Cosine1, sphericalTensor);
makeFunction1Type(Cosine1, symmTensor);
makeFunction1Type(Cosine1, tensor);

makeFunction1Type(Sine1, scalar);
makeFunction1Type(Sine1, vector);
makeFunction1Type(Sine1, sphericalTensor);
makeFunction1Type(Sine1, symmTensor);
makeFunction1Type(Sine1, tensor);

makeFunction1Type(Square1, scalar);
makeFunction1Type(Square1, vector);
makeFunction1Type(Square1, sphericalTensor);
makeFunction1Type(Square1, symmTensor);
makeFunction1Type(Square1, tensor);

} // End namespace Foam


// ************************************************************************* //
