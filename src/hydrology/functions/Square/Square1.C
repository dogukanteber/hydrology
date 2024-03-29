/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "Square1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Square1<Type>::Square1
(
    const word& entryName,
    const dictionary& dict
    #if (OPENFOAM >= 2112)
    , const objectRegistry* obrPtr
    #endif
)
:
    #if (OPENFOAM >= 2112)
    Sine1<Type>(entryName, dict, obrPtr),
    #else
    Sine1<Type>(entryName, dict),
    #endif
    mark_(dict.getOrDefaultCompat<scalar>("mark", {{"markSpace", 2006}}, 1)),
    space_(dict.getOrDefault<scalar>("space", 1))
{}


template<class Type>
Foam::Function1Types::Square1<Type>::Square1(const Square1<Type>& rhs)
:
    Sine1<Type>(rhs),
    mark_(rhs.mark_),
    space_(rhs.space_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Square1<Type>::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("mark", 1, mark_);
    os.writeEntryIfDifferent<scalar>("space", 1, space_);
    Sine1<Type>::writeEntries(os);
}


template<class Type>
void Foam::Function1Types::Square1<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
