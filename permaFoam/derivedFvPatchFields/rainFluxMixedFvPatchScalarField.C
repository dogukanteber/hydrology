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

\*---------------------------------------------------------------------------*/

#include "rainFluxMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::rainFluxMixedFvPatchScalarField::precipitationUnits
>
Foam::rainFluxMixedFvPatchScalarField::unitNames_
{
    { precipitationUnits::METERS_PER_SECOND, "m_s" },
    { precipitationUnits::METERS_PER_HOUR, "m_hour" },
    { precipitationUnits::METERS_PER_MONTH, "m_month" },
    { precipitationUnits::METERS_PER_YEAR, "m_year" },
    { precipitationUnits::MM_PER_SECOND, "mm_s" },
    { precipitationUnits::MM_PER_HOUR, "mm_hour" },
    { precipitationUnits::MM_PER_MONTH, "mm_month" },
    { precipitationUnits::MM_PER_YEAR, "mm_year" },
};


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

#if 0

// The relative hydraulic conductivity, when psi_F < 0

thtil_cl =
(
    pow((1+pow(mag(alpha*psi_F),n)),-(1-(1/n)))
);

Krel_cl =
(
    K
  * pow(thtil_cl,0.5)
  * pow((1-pow((1-pow(thtil_cl,(n/(n-1)))),(1-(1/n)))),2)
);


// When psi_F >= 0, simply use psi = 0 and ignore K, Krel_cl etc
#endif

namespace Foam
{
    // Factor for the relative hydraulic conductivity
    static inline scalar conductivityFactor
    (
        const scalar psi_F,
        const scalar alpha,
        const scalar n
    )
    {
        if (psi_F >= 0)
        {
            return 1;
        }

        const scalar thtil =
            pow((1+pow(mag(alpha*psi_F),n)),-(1-(1/n)));

        return
        (
            sqrt(thtil)
          * sqr(1-pow((1-pow(thtil,(n/(n-1)))),(1-(1/n))))
        );
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static constexpr double sec_per_hour  = (3600.0);
static constexpr double sec_per_day   = (86400.0);
static constexpr double sec_per_year  = (365.0 * sec_per_day);  // common year
static constexpr double sec_per_month = (365.0 * sec_per_day / 12);

namespace Foam
{
    static double getScaleDivisor
    (
        enum rainFluxMixedFvPatchScalarField::precipitationUnits units
    )
    {
        switch (units)
        {
            case rainFluxMixedFvPatchScalarField::METERS_PER_SECOND :
                return 1;

            case rainFluxMixedFvPatchScalarField::METERS_PER_HOUR :
                return sec_per_hour;

            case rainFluxMixedFvPatchScalarField::METERS_PER_MONTH :
                return sec_per_month;

            case rainFluxMixedFvPatchScalarField::METERS_PER_YEAR :
                return sec_per_year;

            case rainFluxMixedFvPatchScalarField::MM_PER_SECOND :
                return 1000;

            case rainFluxMixedFvPatchScalarField::MM_PER_HOUR :
                return 1000 * sec_per_hour;

            case rainFluxMixedFvPatchScalarField::MM_PER_MONTH :
                return 1000 * sec_per_month;

            case rainFluxMixedFvPatchScalarField::MM_PER_YEAR :
                return 1000 * sec_per_year;

            default :
                FatalErrorInFunction
                    << "Programming error - unknown units"
                    << exit(FatalError);
                break;
        }

        return 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rainFluxMixedFvPatchScalarField::
rainFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    units_(precipitationUnits::METERS_PER_SECOND),
    scaleDivisor_(1),
    direction_(0,0,-1),
    facePressureName_(),
    hydConductivityName_(),
    vanGenuchtenName_(),
    invCapillaryLengthName_(),
    rate_(nullptr)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = scalar(1);
}


Foam::rainFluxMixedFvPatchScalarField::
rainFluxMixedFvPatchScalarField
(
    const rainFluxMixedFvPatchScalarField& pfld,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(pfld, p, iF, mapper),
    units_(pfld.units_),
    scaleDivisor_(pfld.scaleDivisor_),
    direction_(pfld.direction_),
    facePressureName_(pfld.facePressureName_),
    hydConductivityName_(pfld.hydConductivityName_),
    vanGenuchtenName_(pfld.vanGenuchtenName_),
    invCapillaryLengthName_(pfld.invCapillaryLengthName_),
    rate_(nullptr)
{}


Foam::rainFluxMixedFvPatchScalarField::
rainFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    units_
    (
        // Use lookupOrDefault instead of getOrDefault for older OpenFOAM
        unitNames_.lookupOrDefault
        (
            "units",
            dict,
            precipitationUnits::METERS_PER_SECOND
        )
    ),
    scaleDivisor_(getScaleDivisor(units_)),
    direction_
    (
        dict.lookupOrDefault<vector>("direction", vector(0,0,-1)).normalise()
    ),
    facePressureName_
    (
        dict.lookupOrDefault<word>("facePressure", "psi_F")
    ),
    hydConductivityName_
    (
        dict.lookupOrDefault<word>("hydConductivity", "K")
    ),
    vanGenuchtenName_
    (
        dict.lookupOrDefault<word>("vanGenuchten", "n")
    ),
    invCapillaryLengthName_
    (
        dict.lookupOrDefault<word>("invCapillaryLength", "alpha")
    ),
    rate_(Function1<scalar>::New("rate", dict))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = scalar(1);
    }
}


Foam::rainFluxMixedFvPatchScalarField::
rainFluxMixedFvPatchScalarField
(
    const rainFluxMixedFvPatchScalarField& pfld,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(pfld, iF),
    units_(pfld.units_),
    scaleDivisor_(pfld.scaleDivisor_),
    direction_(pfld.direction_),
    facePressureName_(pfld.facePressureName_),
    hydConductivityName_(pfld.hydConductivityName_),
    vanGenuchtenName_(pfld.vanGenuchtenName_),
    invCapillaryLengthName_(pfld.invCapillaryLengthName_),
    rate_(pfld.rate_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rainFluxMixedFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vector negativeZ(0,0,-1);

    const label patchi = this->patch().index();
    const polyMesh& mesh = this->patch().boundaryMesh().mesh();

    const auto& sf_psi_F =
        mesh.lookupObject<surfaceScalarField>(facePressureName_);

    const auto& vf_K =
        mesh.lookupObject<volScalarField>(hydConductivityName_);

    const auto& vf_n =
        mesh.lookupObject<volScalarField>(vanGenuchtenName_);

    const auto& vf_alpha =
        mesh.lookupObject<volScalarField>(invCapillaryLengthName_);

    const scalar currRate =
        rate_->value(this->db().time().timeOutputValue()) / scaleDivisor_;


    // Re-zero. Assume fixedValue.
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;

    // Patch fields
    const auto& psi_F = sf_psi_F.boundaryField()[patchi];
    const auto& K = vf_K.boundaryField()[patchi];
    const auto& n = vf_n.boundaryField()[patchi];
    const auto& alpha = vf_alpha.boundaryField()[patchi];

    const auto& Sf = this->patch().Sf();
    const auto& magSf = this->patch().magSf();

    forAll(this->patch(), patchFacei)
    {
        if (psi_F[patchFacei] <= 0)
        {
            const scalar Krel =
            (
                K[patchFacei]
              * conductivityFactor
                (
                    psi_F[patchFacei],
                    alpha[patchFacei],
                    n[patchFacei]
                )
            );

            const scalar gradValue =
            (
                (negativeZ - currRate/Krel * direction_)
              & (Sf[patchFacei] / magSf[patchFacei])
            );

            refGrad()[patchFacei] = gradValue;
            valueFraction()[patchFacei] = 0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::rainFluxMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("units", unitNames_[units_]);
    if (rate_)
    {
        rate_->writeData(os);
    }
    os.writeEntry("direction", direction_);
    os.writeEntry("facePressure", facePressureName_);
    os.writeEntry("hydConductivity", hydConductivityName_);
    os.writeEntry("vanGenuchten", vanGenuchtenName_);
    os.writeEntry("invCapillaryLength", invCapillaryLengthName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    rainFluxMixedFvPatchScalarField
);

} // End namespace Foam


// ************************************************************************* //
