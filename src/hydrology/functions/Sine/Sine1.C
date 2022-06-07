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

#include "Sine1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Sine1<Type>::read(const dictionary& dict)
{
    t0_ = dict.getOrDefault<scalar>("t0", 0);

    clip_.clear();
    dict.readIfPresent("clip", clip_);

    amplitude_.reset(nullptr);
    frequency_.reset(nullptr);
    period_.reset(nullptr);

    // NewIfPresent
    if (dict.found("amplitude", keyType::LITERAL))
    {
        amplitude_ = Function1<scalar>::New("amplitude", dict);
    }
    if (dict.found("period", keyType::LITERAL))
    {
        period_ = Function1<scalar>::New("period", dict);
    }
    if (!period_)
    {
        frequency_ = Function1<scalar>::New("frequency", dict);
    }

    scale_ = Function1<Type>::New("scale", dict);
    level_ = Function1<Type>::New("level", dict);
}


template<class Type>
Foam::Function1Types::Sine1<Type>::Sine1
(
    const word& entryName,
    const dictionary& dict
    #if (OPENFOAM >= 2112)
    , const objectRegistry* obrPtr
    #endif
)
:
    #if (OPENFOAM >= 2112)
    Function1<Type>(entryName, dict, obrPtr),
    #elif (OPENFOAM >= 2012)
    Function1<Type>(entryName, dict),
    #else
    Function1<Type>(entryName),
    #endif
    t0_(0),
    clip_(),
    amplitude_(nullptr),
    frequency_(nullptr),
    period_(nullptr),
    scale_(nullptr),
    level_(nullptr)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::Sine1<Type>::Sine1(const Sine1<Type>& rhs)
:
    Function1<Type>(rhs),
    t0_(rhs.t0_),
    clip_(rhs.clip_),
    amplitude_(rhs.amplitude_.clone()),
    frequency_(rhs.frequency_.clone()),
    period_(rhs.period_.clone()),
    scale_(rhs.scale_.clone()),
    level_(rhs.level_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#if (OPENFOAM >= 2112)
template<class Type>
void Foam::Function1Types::Sine1<Type>::userTimeToTime(const Time& t)
{
    t0_ = t.userTimeToTime(t0_);
}
#endif


template<class Type>
void Foam::Function1Types::Sine1<Type>::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("t0", 0, t0_);
    if (clip_.valid())
    {
        os.writeEntry("clip", clip_);
    }

    if (amplitude_)
    {
        amplitude_->writeData(os);
    }
    if (period_)
    {
        period_->writeData(os);
    }
    if (frequency_)
    {
        frequency_->writeData(os);
    }
    scale_->writeData(os);
    level_->writeData(os);
}


template<class Type>
void Foam::Function1Types::Sine1<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
