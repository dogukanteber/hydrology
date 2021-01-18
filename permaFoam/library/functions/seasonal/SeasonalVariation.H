/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::Function1Types::SeasonalVariation

Description

    Templated function for mimicking seasonal meteorological variations
    as the combination of functions representing two variation modes
    for a given variable.

    \code
        pos0(selector) * mode1 + neg(selector) * mode2
    \code

    The scalar selector function is used to determine which mode is
    active at a particular time.

    The dictionary specification could resemble this:
    \verbatim
    entry1
    {
        type        seasonal;
        dualMode    true;  // default

        modeSelector
        {
           type        sine;
           t0          27594000.0;  // 319 days
           period      31536000.0;  // 365 days
           scale       1;
           level       0.2
        }

        mode1
        {
            type        sine1;
            t0          22170240.0;  // 256.6 days
            period      31536000.0;  // 365 days
            scale       20.3;
            level       265;
            clip        (258.43 1000);  // min -15C
        }

        mode2
        {
            type        sine1;
            t0          15768000.0;  // 182.5 days
            period      31536000.0;  // 365 days

            scale       13.365;
            level       267.31;
            clip        (258.43 1000);  // min -15C
        }
    }
    \endverbatim

    where the entries mean:
    \table
        Property  | Description             | Type | Reqd | Default
        type      | Function type: seasonal | word | yes |
        dualMode  | Use mode2 as well as mode1 | Switch | no | ON
        modeSelector | Select mode1 or mode2 | Function1<scalar> | conditional |
        mode1     | Values for mode1        | Function1<Type> | yes |
        mode2     | Values for mode2        | Function1<Type> | conditional |
    \endtable

SourceFiles
    SeasonalVariationI.H
    SeasonalVariation.C

\*---------------------------------------------------------------------------*/

#ifndef function1Types_SeasonalVariation_H
#define function1Types_SeasonalVariation_H

#include "Function1.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                      Class SeasonalVariation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SeasonalVariation
:
    public Function1<Type>
{
    // Private Data

        //- Enable use of dual mode (default: on)
        Switch dualMode_;

        //- Selector for mode1/mode2
        autoPtr<Function1<scalar>> selector_;

        //- Mode 1
        autoPtr<Function1<Type>> mode1_;

        //- Mode 2
        autoPtr<Function1<Type>> mode2_;


    // Private Member Functions

        //- Read coefficients dictionary
        void read(const dictionary& dict);

public:

    // Runtime type information
    TypeName("seasonal");


    // Generated Methods

        //- No copy assignment
        void operator=(const SeasonalVariation<Type>&) = delete;


    // Constructors

        //- Construct from entry name and dictionary
        SeasonalVariation(const word& entryName,const dictionary& dict);

        //- Copy construct
        explicit SeasonalVariation(const SeasonalVariation<Type>& rhs);


    //- Destructor
    virtual ~SeasonalVariation() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SeasonalVariationI.H"

#ifdef NoRepository
    #include "SeasonalVariation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
