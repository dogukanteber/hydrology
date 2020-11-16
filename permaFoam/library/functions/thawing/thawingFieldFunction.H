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

Class
    Foam::PatchFunction1Types::thawingFieldFunction

Description
    A scalar PatchFunction1 to decide if patch faces are thawing.

    \f[
        (Tface - (Tmelt + Tshift)) >= 0
    \f]

    where

    \vartable
        Name       | Description                      | Data type
        type       | thawingField                     |
        Tface      | Temperature field on faces [K]   | surfScalarField
        Tmelt      | Temperature for phase change [K] | constant/volScalarField
        Tshift     | thermal shift [K]                | constant (optional)
    \endvartable

    Note the \c Tmelt can be specified as a constant value or a field name.

SourceFiles
    thawingFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef function1s_thawingFieldFunction_H
#define function1s_thawingFieldFunction_H

#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                    Class thawingFieldFunction Declaration
\*---------------------------------------------------------------------------*/

class thawingFieldFunction
:
    public PatchFunction1<scalar>
{
    // Private Data

        //- Tface : temperature face field name
        word faceName_;

        //- Tmelt : melting temperature volume field (optional)
        word meltName_;

        //- Tmelt : melting temperature value
        scalar meltValue_;

        //- Tshift : shift for stabilization
        scalar shiftValue_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const thawingFieldFunction&) = delete;


public:

    // Runtime type information
    TypeName("thawingField");


    // Constructors

        //- Construct from entry name and dictionary
        thawingFieldFunction
        (
            const polyPatch& pp,
            const word& type,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct
        explicit thawingFieldFunction(const thawingFieldFunction& rhs);

        //- Copy construct, setting patch
        explicit thawingFieldFunction
        (
            const thawingFieldFunction& rhs,
            const polyPatch& pp
        );

        //- Construct and return a clone
        virtual tmp<PatchFunction1<scalar>> clone() const
        {
            return tmp<PatchFunction1<scalar>>
            (
                new thawingFieldFunction(*this)
            );
        }

        //- Construct and return a clone setting patch
        virtual tmp<PatchFunction1<scalar>> clone(const polyPatch& pp) const
        {
            return tmp<PatchFunction1<scalar>>
            (
                new thawingFieldFunction(*this, pp)
            );
        }


    //- Destructor
    virtual ~thawingFieldFunction() = default;


    // Member Functions

    // Evaluation

        //- Return a field of 0/1 as false/true or thaw factor
        virtual tmp<Field<scalar>> value(const scalar x) const;

        //- Is value constant (i.e. independent of x)
        virtual inline bool constant() const
        {
            return false;
        }

        //- Is value uniform (i.e. independent of coordinate)
        virtual inline bool uniform() const
        {
            return false;
        }


    // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const FieldMapper& mapper);

        //- Reverse map the given PatchFunction1 onto this PatchFunction1
        virtual void rmap
        (
            const PatchFunction1<scalar>& pf1,
            const labelList& addr
        );


    // IO

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
