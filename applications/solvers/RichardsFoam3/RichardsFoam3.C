/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 Laurent Orgogozo
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

Application
    RichardsFoam3

Description
    Transient solver for water transfers in unsaturated porous
        media.
    See Orgogozo et al. CPC 2014 and Orgogozo CPC 2015 for the implementation
        and validation of water transfers resolution.
    Handling of non-linearities : Picard loops.
    Coupling between equations : Sequential operator splitting.
    Global computations of the convergence criteria.
    Empirical adaptative time stepping with stabilisation procedure.
    NB: use backward scheme for time discretisation.

Contributors to the design and the writing of the code (chronological order):
    Laurent Orgogozo
    Cyprien Soulaine
    Michel Quintard
    Rachid Ababou
    Mark Olesen

References
    [1] Orgogozo L. 2015. RichardsFoam2: a new version of RichardsFoam devoted
      to the modelling of the vadose zone. Computer Physics Communications,
      196, 619-620.
      https://doi.org/10.1016/j.cpc.2015.07.009
      https://hal.archives-ouvertes.fr/hal-01299854
    [2] Orgogozo L., Renon N., Soulaine C., Hénon F., Tomer S.K., Labat D.,
      Pokrovsky O.S., Sekhar M., Ababou R., Quintard M. 2014a. An open source
      massively parallel solver for Richards equation: Mechanistic modelling of
      water fluxes at the watershed scale. Computer Physics Communications,
      185, 3358-3371.
      https://doi.org/10.1016/j.cpc.2014.08.004
      http://oatao.univ-toulouse.fr/12140/

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "declarePicardControls.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//                           Initialisation
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    // Number of mesh cells

    const label nbMesh = returnReduce(mesh.nCells(), sumOp<label>());

    // Initialisation of the Richards equation parameters

    const volScalarField z(mesh.C().component(vector::Z));

    volScalarField psi_tmp = psi;

    volScalarField psim1 = psi;

    // Equation (2), in RichardsFoam3_systemOfEquations.pdf
    volScalarField thtil
    (

                           pos0(psi)
                           +
                           neg(psi)*
                           pow(
                                  (1 + pow(
                                              mag(alpha*psi),
                                              n
                                          )
                                   ),
                                  - (1 - (1/n)))
    );


    // Equation (2), in RichardsFoam3_systemOfEquations.pdf
    volScalarField thtil_tmp
    (
                               pos0(psi_tmp)+neg(psi_tmp)*
                               pow(
                                      (1 + pow(
                                                  mag(alpha*psi_tmp),
                                                  n
                                              )
                                      ),
                                      -(1 - (1/n)))
    );


    // Equation (4), in RichardsFoam3_systemOfEquations.pdf
    volScalarField Krel
    (
                          (
                                      pos0(psi)*
                                      K
                                      +
                                      neg(psi)*
                                      K*
                                      pow(thtil,0.5)*
                                      pow(
                                             (1-pow(
                                                      (1-pow(thtil,(n/(n-1)))),
                                                      (1 - (1/n)))
                                                   ),
                                             2
                                         )
                          )
    );

    // Equation (3), in RichardsFoam3_systemOfEquations.pdf
    volScalarField Crel
    (
                          S +     neg(psi)*
                                  (
                                      (thetas - thetar)*
                                      (thtil - thtil_tmp)*
                                      (
                                          1./
                                          (
                                              (
                                                  usf*
                                                  pos0(psi - psi_tmp)*
                                                  pos0(psi_tmp - psi)
                                              )
                                              +
                                              psi
                                              -
                                              psi_tmp
                                          )
                                      )
                                  )
    );

    volVectorField gradk(fvc::grad(Krel));

    volScalarField gradkz(gradk.component(vector::Z));

// Initialisation of the velocity field computation

    surfaceScalarField flupsi(fvc::snGrad(psi));

    surfaceScalarField fluK(fvc::interpolate(Krel,"interpolate(Krel)"));

    surfaceVectorField norFace(mesh.Sf()/mesh.magSf());

    surfaceScalarField fluvuz(norFace.component(vector::Z));

    U = -fluK*(flupsi + fluvuz);

    Uvol =
      (
          -Krel * (fvc::grad(psi,"grad(psi)") + vuz)
      );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//                    Beginning of the resolution
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        // Loop 1 opening - temporal iterations //
        #include "readPicardControls.H"
        #include "updatePicardControls.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Residual tests for exiting the Picard loops
        scalar convergeFlow = 0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//                    Richards equation solving
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        for (int cyc = 0; cyc < nMaxCycle; ++cyc)
        {
            // Loop 2 opening - adaptive time stepping loop //

// Computation of the Actual Evapotranspiration sink term
// for the current time step

            // Equation (5), in RichardsFoam3_systemOfEquations.pdf
            {

            volScalarField quotient
                 (
                     (theta - thetaWP)
                   / (runTime.deltaTValue()*tdim)
                 );

            AET = PET*
                  pos0(quotient - PET)
                  +
                  quotient*
                  neg(quotient - PET)*
                  neg(thetaWP-theta);
            }

            for (int picFlowi = 0; picFlowi < nIterPicardFlow; ++picFlowi)
            {
                // Loop 3 opening - Picard iterations for Richards equation //

                psim1 = psi;

                psi = psi_tmp;

// Resolution of the linear system.

                // Equation (1), in RichardsFoam3_systemOfEquations.pdf
                {
                    fvScalarMatrix psiEqn
                    (
                        Crel*fvm::ddt(psi)
                     == fvm::laplacian(Krel, psi, "laplacian(Krel,psi)")
                      + gradkz
                      - AET
                    );
                    psiEqn.solve();
                }

// update of the varying transport properties.

                // Equation (2), in RichardsFoam3_systemOfEquations.pdf
                thtil =     pos0(psi) + neg(psi)*
                            pow(
                                   (1 + pow(
                                                mag(alpha*psi),
                                                n
                                            )),
                                   -(1 - (1/n))
                         );


                // Equation (4), in RichardsFoam3_systemOfEquations.pdf
                Krel = (
                               pos0(psi)*K
                               +
                               neg(psi)*
                               K*
                               pow(thtil, 0.5)*
                               pow(
                                      (1 - pow(
                                                  (1 -pow(thtil, (n/(n-1)))),
                                                  (1-(1/n))
                                              )),
                                      2
                                  )
                       );

                // Equation (3), in RichardsFoam3_systemOfEquations.pdf
                Crel = S +         neg(psi)*
                                   (
                                       (thetas - thetar)*
                                       (thtil - thtil_tmp)*
                                       (1./((usf*pos0(psi - psi_tmp)*pos0(
                                              psi_tmp - psi)) + psi - psi_tmp))
                                   );

// update of the gravity term.
                gradk = fvc::grad(Krel);

                gradkz = gradk.component(vector::Z);

// updated of water content
                theta = (thetas - thetar)*thtil + thetar;

// exit test of the Picard loop.
// note the use of gSum instead of sum.
                err = psi - psim1;

                convergeFlow = gSumMag(err)/nbMesh;

                if (convergeFlow < precPicardFlow)
                {
                    currentPicardFlow = picFlowi;
                    Info<< "Error psi = " << convergeFlow
                        << " Picard Flow = " << currentPicardFlow
                        << nl << endl;
                    break;
                }

            } // loop 3 closing - Picard iterations for Richards equation //


// Update of the velocity field
            flupsi = fvc::snGrad(psi);

            fluK = fvc::interpolate(Krel,"interpolate(Krel)");

            U = - fluK*(flupsi + fluvuz);

            Uvol =
            (
              - Krel * (fvc::grad(psi,"grad(psi)") + vuz)
            );


// exit test of the temporal iteration loop ; exit if Richards
// resolution is succesfull

            if
            (
                (convergeFlow < precPicardFlow)
            )
            {
                // Converged
                break;
            }

            bool needTimeAdjustment = false;
            if
            (
                (convergeFlow >= precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Richards criterium not reached"
                    " - computation restarted with a smaller time step. "
                    << "Error psi = " << convergeFlow << nl;
            }

            if (needTimeAdjustment)
            {
                runTime.setDeltaT((1/tFact)*runTime.deltaTValue());
                Info<< "deltaT = " <<  runTime.deltaTValue() << nl << endl;
            }

        } // Loop 2 closing - adaptive time stepping loop //


        // Display error messages, if any
        {
            bool needTimeAdjustment = false;
            if
            (
                (convergeFlow >= precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Failure in Richards convergence / Error psi = "
                    << convergeFlow << nl;
            }

            if (needTimeAdjustment)
            {
                runTime.setDeltaT((1/tFact)*runTime.deltaTValue());
                Info<< "deltaT = " <<  runTime.deltaTValue() << nl << endl;
            }
        }


// Storage of the field of results for the Picard loops of the next time step
        psi_tmp = psi;
        psi_F = fvc::interpolate(psi,"interpolate(psi)");
        H = psi + z;
        // Equation (2), in RichardsFoam3_systemOfEquations.pdf
        thtil_tmp =     pos0(psi)
                        +
                        neg(psi)*
                        pow(
                               (1 + pow(mag(alpha*psi_tmp),n)),
                               - (1 - (1/n))
                    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//                    Writing of the output results
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


        runTime.write();

        runTime.printExecutionTime(Info);

    } // Loop 1 closing - temporal loop //

    Info<< "\nEnd\n" << endl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
