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

#include "SeasonalVariation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::SeasonalVariation<Type>::read
(
    const dictionary& dict
)
{
    dualMode_ = dict.getOrDefault("dualMode", Switch(Switch::ON));

    selector_.reset(nullptr);
    mode2_.reset(nullptr);

    mode1_ = Function1<Type>::New("mode1", dict);

    if (dualMode_)
    {
        selector_ = Function1<scalar>::New("modeSelector", dict);
        mode2_ = Function1<Type>::New("mode2", dict);
    }
    else
    {
        // NewIfPresent

        if (dict.found("modeSelector", keyType::LITERAL))
        {
            selector_ = Function1<scalar>::New("modeSelector", dict);
        }
        if (dict.found("mode2", keyType::LITERAL))
        {
            mode2_ = Function1<Type>::New("mode2", dict);
        }
    }
}


template<class Type>
Foam::Function1Types::SeasonalVariation<Type>::SeasonalVariation
(
    const word& entryName,
    const dictionary& dict
)
:
    #if (OPENFOAM >= 2012)
    Function1<Type>(entryName, dict),
    #else
    Function1<Type>(entryName),
    #endif
    dualMode_(Switch::ON),
    selector_(nullptr),
    mode1_(nullptr),
    mode2_(nullptr)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::SeasonalVariation<Type>::SeasonalVariation
(
    const SeasonalVariation<Type>& rhs
)
:
    Function1<Type>(rhs),
    dualMode_(rhs.dualMode_),
    selector_(rhs.selector_.clone()),
    mode1_(rhs.mode1_.clone()),
    mode2_(rhs.mode2_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::SeasonalVariation<Type>::writeData
(
    Ostream& os
) const
{
    if (dualMode_)
    {
        os.writeEntry("dualMode", dualMode_);
    }

    if (selector_)
    {
        selector_->writeData(os);
    }

    mode1_->writeData(os);

    if (mode2_)
    {
        mode2_->writeData(os);
    }
}


// ************************************************************************* //
