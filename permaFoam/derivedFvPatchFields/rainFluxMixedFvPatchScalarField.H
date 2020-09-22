/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd
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
    Foam::rainFluxMixedFvPatchScalarField

Description
    Mixed boundary condition for rainfall.

Usage
    \table
        Property     | Description                  | Required | Default
        type         | "rainFallFlux"               | yes |
        rate         | Rainfall rate                | yes |
        units        | Rainfall units               | no  | m_s
        direction    | Direction of rain            | no  | (0 0 -1)
        facePressure | Name of face pressure field  | no  | psi_F
        hydConductivity | Name of hydraulic conductivity | no  | K
        vanGenuchten | Name of van Genuchten        | no  | n
        invCapillaryLength | Name of inverse capillary length | no  | alpha
    \endtable

    The rainfall units
    \table
        Name        | Description
        m_s         | Meters/second
        m_hour      | Meters/hour
        m_month     | Meters/month
        m_year      | Meters/year
        mm_s        | Millimeters/second
        mm_hour     | Millimeters/hour
        mm_month    | Millimeters/month
        mm_year     | Millimeters/year
    \endtable

SourceFiles
    rainFluxMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rainFluxMixedFvPatchScalarField_H
#define rainFluxMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Enum.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class rainFluxMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rainFluxMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

    // Enumerations

        //- Units for precipitation
        enum precipitationUnits
        {
            // Meters per time
            METERS_PER_SECOND = 1,
            METERS_PER_HOUR,
            METERS_PER_MONTH,
            METERS_PER_YEAR,

            // Millimeters per time
            MM_PER_SECOND,
            MM_PER_HOUR,
            MM_PER_MONTH,
            MM_PER_YEAR
        };

private:

    // Private Data

        //- Names for precipitation units
        static const Enum<precipitationUnits> unitNames_;

        //- Precipitation units
        const precipitationUnits units_;

        //- Scale divisor for precipitation units
        const scalar scaleDivisor_;

        //- Direction of rainfall: default (0 0 -1)
        const vector direction_;

        //- Name of water pressure at faces field
        const word facePressureName_;

        //- Name of saturated hydraulic conductivity field
        const word hydConductivityName_;

        //- Name of van Genuchten exponent field
        const word vanGenuchtenName_;

        //- Name of inverse of capillary length field
        const word invCapillaryLengthName_;

        //- The rainfall rate
        autoPtr<Function1<scalar>> rate_;


public:

    //- Runtime type information
    TypeName("rainFallFlux");


    // Constructors

        //- Construct from patch and internal field
        rainFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rainFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rainFluxMixedFvPatchScalarField
        (
            const rainFluxMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new rainFluxMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        rainFluxMixedFvPatchScalarField
        (
            const rainFluxMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new rainFluxMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
