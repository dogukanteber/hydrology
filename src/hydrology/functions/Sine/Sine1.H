/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::Function1Types::Sine1

Description
    A templated sine function, with support for offset, clipping etc.

    The wave period can be specified directly

        \f[
            clip(a sin(2 \pi (t - t_0) / p)) s + l)
        \f]

    Or it can be specified by the frequency

        \f[
           clip(a sin(2 \pi f (t - t_0)) s + l)
        \f]

    where
    \vartable
        Symbol  | Description       | Units
        a       | Amplitude         | -
        f       | Frequency         | [1/s]
        p       | Period            | [s]
        s       | Type scale factor | -
        l       | Type offset level | -
        clip    | clipping range    | -
        t       | Time              | [s]
        t_0     | Start time offset | [s]
    \endvartable

    The dictionary specification would typically resemble this:
    \verbatim
    entry1
    {
        type        sine1;
        frequency   10;
        amplitude   0.1;

        // A scalar Function1
        scale       2e-6;
        level       2e-6;
    }
    entry2
    {
        type        sine1;
        frequency   10;

        // A vector Function1
        scale       (1 0.1 0);
        level       (10 1 0);
    }
    \endverbatim

    where the entries mean:
    \table
        Property  | Description             | Type | Reqd | Default
        type      | Function type: sine1    | word | yes |
        amplitude | Amplitude               | Function1<scalar> | no | 1
        frequency | Frequency [1/s]         | Function1<scalar> | or period |
        period    | Period [s]              | Function1<scalar> | or frequency |
        scale     | Scale factor (Type)     | Function1<Type> | yes |
        level     | Offset level (Type)     | Function1<Type> | yes |
        clip      | Clipping (Type)         | MinMax<Type> | no | none
        t0        | Start time offset       | scalar       | no | 0
    \endtable

Note
    For slow oscillations it can be more intuitive to specify the period.

SourceFiles
    Sine1.C
    Sine1I.H

\*---------------------------------------------------------------------------*/

#ifndef function1Types_Sine1_H
#define function1Types_Sine1_H

#include "Function1.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                            Class Sine Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Sine1
:
    public Function1<Type>
{
protected:

    // Protected Data

        //- Start-time for the function
        scalar t0_;

        //- Clipping to apply to the output
        MinMax<Type> clip_;

        //- Scalar amplitude of the function (optional)
        autoPtr<Function1<scalar>> amplitude_;

        //- Frequency of the function (or specify period)
        autoPtr<Function1<scalar>> frequency_;

        //- Period of the function (or specify frequency)
        autoPtr<Function1<scalar>> period_;

        //- Scaling factor for the function
        autoPtr<Function1<Type>> scale_;

        //- Level to add to the scaled function
        autoPtr<Function1<Type>> level_;


    // Protected Member Functions

        //- Read coefficients from the dictionary
        void read(const dictionary& dict);

        //- The cycle: (freq * time) or (time / period)
        inline scalar cycle(const scalar t) const;

        //- Calculated cos value at time t
        inline scalar cosForm(const scalar t) const;

        //- Calculated sin value at time t
        inline scalar sinForm(const scalar t) const;

        //- Calculated square value at time t.
        //  The positive fraction is 0-1
        inline scalar squareForm(const scalar t, const scalar posFrac) const;

        //- Return value for time t, using cos form
        inline Type cosValue(const scalar t) const;

        //- Return value for time t, using sin form
        inline Type sinValue(const scalar t) const;

        //- Return value for time t, using square form
        inline Type squareValue(const scalar t, const scalar posFrac) const;


public:

    // Runtime type information
    TypeName("sine1");


    // Generated Methods

        //- No copy assignment
        void operator=(const Sine1<Type>&) = delete;


    // Constructors

        //- Construct from entry name and dictionary
        Sine1(const word& entryName, const dictionary& dict);

        //- Copy construct
        explicit Sine1(const Sine1<Type>& rhs);


    //- Destructor
    virtual ~Sine1() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const
        {
            return Sine1<Type>::sinValue(t);
        }

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Sine1I.H"

#ifdef NoRepository
    #include "Sine1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
