/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::Function1Types::Cosine

Description
    A templated cosine function, with support for offset etc.

    The wave period can be specified directly

        \f[
            clip(a cos(2 \pi (t - t_0) / p)) s + l)
        \f]

    Or it can be specified by the frequency

        \f[
           clip(a cos(2 \pi f (t - t_0)) s + l)
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
        t0      | Start time offset | [s]
    \endvartable

    The dictionary specification would typically resemble this:
    \verbatim
    entry1
    {
        type        cosine;
        frequency   10;
        amplitude   0.1;

        // A scalar Function1
        scale       2e-6;
        level       2e-6;
    }
    entry2
    {
        type        cosine;
        frequency   10;

        // A vector Function1
        scale       (1 0.1 0);
        level       (10 1 0);
    }
    \endverbatim

    where the entries mean:
    \table
        Property  | Description             | Type | Reqd | Default
        type      | Function type: cosine1  | word | yes |
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

\*---------------------------------------------------------------------------*/

#ifndef hydrology_function1Types_Cosine1_H
#define hydrology_function1Types_Cosine1_H

#include "Sine1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                            Class Cosine Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Cosine1
:
    public Function1Types::Sine1<Type>
{
public:

    // Runtime type information
    TypeName("cosine1");


    // Generated Methods

        //- No copy assignment
        void operator=(const Cosine1<Type>&) = delete;


    // Constructors

        //- Construct from entry name and dictionary
        Cosine1
        (
            const word& entryName,
            const dictionary& dict
            #if (OPENFOAM >= 2112)
            , const objectRegistry* obrPtr = nullptr
            #endif
        )
        :
            #if (OPENFOAM >= 2112)
            Sine1<Type>(entryName, dict, obrPtr)
            #else
            Sine1<Type>(entryName, dict)
            #endif
        {}

        //- Copy construct
        explicit Cosine1(const Cosine1<Type>& rhs)
        :
            Sine1<Type>(rhs)
        {}

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Cosine1<Type>(*this));
        }


    //- Destructor
    virtual ~Cosine1() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const
        {
            return Sine1<Type>::cosValue(t);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
