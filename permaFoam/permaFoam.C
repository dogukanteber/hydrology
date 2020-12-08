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
    permaFoam

Description
    Transient solver for water and thermal transfers in unsaturated porous
        media with freeze/thaw of the pore water.
    See Orgogozo et al. PPP 2019 for the implementation of thermal transfers
        resolution.
    See Grenier et al. AWR 2018 for the validation of the numerical resolution
        thermal transfers.
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
    Christophe Grenier
    Mark Olesen

References
    [1] Orgogozo, L., Prokushkin, A.S., Pokrovsky, O.S., Grenier, C., Quintard,
      M., Viers, J., Audry, S. (2019). Water and energy transfer modeling in a
      permafrost-dominated, forested catchment of Central Siberia: the key role
      of rooting depth. Permafrost and Periglacial Processes, 30, 75-89.
      https://doi.org/10.1002/ppp.1995
      https://hal.archives-ouvertes.fr/hal-02014619
    [2] Grenier, C., Anbergen, H., Bense, V., Chanzy, Q., Coon, E., Collier,
      N., Costard, F., Ferry, M., Frampton, A., Frederick, J., Gonçalvès, J.,
      Holmén, J., Jost, A., Kokh, S., Kurylyk, B., McKenzie, J., Molson, J.,
      Mouche, E., Orgogozo, L., Pannetier, R., Rivière, A., Roux, N., Rühaak,
      W., Scheidegger, J., Selroos J.-O., Therrien R., Vidstrand P., Voss C.
      (2018). Groundwater flow and heat transport for systems undergoing
      freeze-thaw: Intercomparison of numerical simulators for 2D test cases.
      Advances in Water Resources, 114, 196-218.
      https://doi.org/10.1016/j.advwatres.2018.02.001
      https://hal.archives-ouvertes.fr/hal-01396632
    [3] Orgogozo L. 2015. RichardsFoam2: a new version of RichardsFoam devoted
      to the modelling of the vadose zone. Computer Physics Communications,
      196, 619-620.
      https://doi.org/10.1016/j.cpc.2015.07.009
      https://hal.archives-ouvertes.fr/hal-01299854
    [4] Orgogozo L., Renon N., Soulaine C., Hénon F., Tomer S.K., Labat D.,
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

// Initialisaton of the test of non-advection of ice

    volScalarField fTestConv(theta - thetag);


    // Initialisation of the Richards equation parameters

    const volScalarField z(mesh.C().component(vector::Z));

    volScalarField psi_tmp = psi;

    volScalarField psim1 = psi;

    // Equation (A.5), in [1] Appendix S1
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


    // Equation (A.5), in [1] Appendix S1
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


    // Equations (A.4), (A.7) and (A.8), in [1] Appendix S1
    volScalarField Krel
    (
                          max(KrThMin,pow(10.,-omega*thetag))*
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

    // Equation (A.6), in [1] Appendix S1
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

    // Equation (A.2), in [1] Appendix S1
    U = -fluK*(flupsi + fluvuz);

    surfaceScalarField phi(Cthw*U*mesh.magSf());

    // Equation (A.2), in [1] Appendix S1
    Uvol =
      (
          -Krel * (fvc::grad(psi,"grad(psi)") + vuz)
      );

// Initialisation of the Heat transfer equation parameters

    // Equations (A.3) and (A.9), in [1] Appendix S1
    volScalarField Cdg
    (
                       - L*
                         pos0(Tmelt - T)*
                         (theta - thetar)*
                         (2.*(T - Tmelt)/(W*W))*
                         exp(
                                -pow(((T - Tmelt)/W),2)
                            )
    );

    volScalarField T_tmp = T;

    volScalarField Tm1 = T;

    // Equation (A.13), in [1] Appendix S1
    volScalarField Kth
    (
         Kthdim
       * pow(Kthsoil*(1./Kthdim),(1. - thetas))
       * pow(Kthw*(1./Kthdim),thetal)
       * pow(Kthice*(1./Kthdim),thetag)
       * pow(Kthair*(1./Kthdim),(thetas - theta))
    );

    // Equation (A.14), in [1] Appendix S1
    volScalarField Cth
    (
          (1. - thetas)*Cthsoil
        + thetal*Cthw
        + thetag*Cthice
        + (thetas - theta)*Cthair
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
        scalar convergeThermal = 0;
        scalar testConv = 0;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//                    Richards equation solving
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        for (int cyc = 0; cyc < nMaxCycle; ++cyc)
        {
            // Loop 2 opening - adaptive time stepping loop //

// Computation of the Actual Evapotranspiration sink term
// for the current time step

            // Equations (A.11) and (A.12), in [1] Appendix S1
            {

            volScalarField quotient
                 (
                     (thetal - thetaWP)
                   / (runTime.deltaTValue()*tdim)
                 );

            AET = PET*
                  pos0(quotient - PET)
                  +
                  quotient*
                  neg(quotient - PET)*
                  neg(thetaWP-thetal);
            }

            for (int picFlowi = 0; picFlowi < nIterPicardFlow; ++picFlowi)
            {
                // Loop 3 opening - Picard iterations for Richards equation //

                psim1 = psi;

                psi = psi_tmp;

// Resolution of the linear system.

                // Equation (A.1), in [1] Appendix S1
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

                // Equation (A.5), in [1] Appendix S1
                thtil =     pos0(psi) + neg(psi)*
                            pow(
                                   (1 + pow(
                                                mag(alpha*psi),
                                                n
                                            )),
                                   -(1 - (1/n))
                         );


                // Equations (A.4), (A.7) and (A.8), in [1] Appendix S1
                Krel = max(KrThMin, pow(10., - omega*thetag))*
                       (
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

                // Equation (A.6), in [1] Appendix S1
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

// test of non-advection of ice
                fTestConv = neg(theta - thetag);

                testConv = gSum(fTestConv);

// exit test of the Picard loop.
// note the use of gSum instead of sum.
                err = psi - psim1;

                convergeFlow = gSumMag(err)/nbMesh;

                if ((convergeFlow < precPicardFlow) && equal(testConv, 0))
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

            // Equation (A.2), in [1] Appendix S1
            U = - fluK*(flupsi + fluvuz);

            // Equation (A.1), in [1] Appendix S1
            Uvol =
            (
              - Krel * (fvc::grad(psi,"grad(psi)") + vuz)
            );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//                    Heat transfer equation solving
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// thermal fluxes normal to faces
            phi = U*mesh.magSf()*Cthw;

// Update of the fields of apparent thermal properties
            // Equation (A.13), in [1] Appendix S1
            Kth =
            (
                Kthdim
              * pow(Kthsoil*(1./Kthdim),(1. - thetas))
              * pow(Kthw*(1./Kthdim),thetal)
              * pow(Kthice*(1./Kthdim),thetag)
              * pow(Kthair*(1./Kthdim),(thetas - theta))
            );

            // Equation (A.14), in [1] Appendix S1
            Cth =
            (
                (1. - thetas)*Cthsoil
              + thetal*Cthw
              + thetag*Cthice
              + (thetas - theta)*Cthair
            );

            for (int picThermi = 0; picThermi < nIterPicardThermal; ++picThermi)
            {
                // loop 4 opening - Picard iteration for thermal equation //

                Tm1 = T;

                T = T_tmp;

//  Resolution of thermal transfers
                // Equation (A.3), in [1] Appendix S1
                {
                    fvScalarMatrix thEqn
                    (
                        fvm::ddt(Cth+Cdg, T)
                      + fvm::div(phi, T, "div(phi,T)")
                      - fvm::laplacian(Kth, T, "laplacian(Kth,T)")
                    );
                    thEqn.solve();
                }

// Application of forced temperature bounds
                T = min(max(T, Tminc), Tmaxc);

//  Computation of phase changes
                // Equation (A.9), in [1] Appendix S1
                thetag = max(
                                0.,
                                pos0(Tmelt - T)*
                                    (
                                        theta
                                        -
                                        (
                                            (theta - thetar)*exp( - pow(((T -
                                                                 Tmelt)/W), 2))
                                          + thetar
                                        )
                                    )
                            );

                // Equation (A.10), in [1] Appendix S1
                thetal = theta - thetag;

// Update of the fields of apparent thermal properties

                // Equation (A.13), in [1] Appendix S1
                Kth =
                (
                    Kthdim
                  * pow(Kthsoil/Kthdim, (1. - thetas))
                  * pow(Kthw/Kthdim, thetal)
                  * pow(Kthice/Kthdim, thetag)
                  * pow(Kthair/Kthdim, (thetas - theta))
                );

                // Equation (A.14), in [1]  Appendix S1
                Cth =
                (
                    (1. - thetas)*Cthsoil
                  + thetal*Cthw
                  + thetag*Cthice
                  + (thetas - theta)*Cthair
                );


// exit test of the thermal Picard loop
// note the use of gSum instead of sum.

                errT = mag(T - Tm1);

                convergeThermal = gMax(errT);

                if (convergeThermal < precPicardThermal)
                {
                    currentPicardThermal = picThermi;
                    Info<< "Error T = " << convergeThermal
                        << " Picard Thermal = " << currentPicardThermal
                        << nl << endl;
                    break;
                }

            } // loop 4 closing - Picard iterations for thermal equation //


// exit test of the temporal iteration loop ; exit if Richards AND thermal
// resolutions are succesfull, and if there is no detectable ice advection

            if
            (
                (convergeFlow < precPicardFlow)
             && (convergeThermal < precPicardThermal)
             && equal(testConv, 0)
            )
            {
                // Converged
                break;
            }

            bool needTimeAdjustment = false;
            if
            (
                (convergeThermal >= precPicardThermal)
             && (convergeFlow < precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Thermal criterium not reached"
                    " - computation restarted with a smaller time step. "
                    << "Error T = " << convergeThermal << nl;
            }
            else if
            (
                (convergeThermal < precPicardThermal)
             && (convergeFlow >= precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Richards criterium not reached"
                    " - computation restarted with a smaller time step. "
                    << "Error psi = " << convergeFlow << nl;
            }
            else if (testConv > 0.)
            {
                needTimeAdjustment = true;
                Info<< "Ice advection"
                    " - computation restarted with a smaller time step. "
                    << "Affected cells = " << testConv << nl;
            }
            else
            {
                needTimeAdjustment = true;
                Info<< "Flow and Thermal criteria not reached"
                    " - computation restarted with a smaller time step. "
                    << nl
                    << "Error psi = " << convergeFlow << nl
                    << "Error T = " << convergeThermal << nl;
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
                (convergeThermal >= precPicardThermal)
             && (convergeFlow < precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Failure in thermal convergence / Error T = "
                    << convergeThermal << nl;
            }
            else if
            (
                (convergeThermal < precPicardThermal)
             && (convergeFlow >= precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Failure in Richards convergence / Error psi = "
                    << convergeFlow << nl;
            }
            else if
            (
                (convergeThermal >= precPicardThermal) &&
                (convergeFlow >= precPicardFlow)
            )
            {
                needTimeAdjustment = true;
                Info<< "Failure in Richards and thermal convergence" << nl
                    << "Error psi = " << convergeFlow << nl
                    << "Error T = " << convergeThermal << nl;
            }

            if (needTimeAdjustment)
            {
                runTime.setDeltaT((1/tFact)*runTime.deltaTValue());
                Info<< "deltaT = " <<  runTime.deltaTValue() << nl << endl;
            }
        }


//  Update of the latent heat part of the apparent thermal capacity
        // Equations (A.3) and (A.9), in [1] Appendix S1
        Cdg = -
              L*
              pos0(Tmelt - T)*
              (theta-thetar)*
              (2.*(T - Tmelt)/(W*W))*
              exp(-pow(((T - Tmelt)/W),2));

// Storage of the field of results for the Picard loops of the next time step
        psi_tmp = psi;
        psi_F = fvc::interpolate(psi,"interpolate(psi)");
        // Equation (A.5), in [1] Appendix S1
        thtil_tmp =     pos0(psi)
                        +
                        neg(psi)*
                        pow(
                               (1 + pow(mag(alpha*psi_tmp),n)),
                               - (1 - (1/n))
                    );
        T_tmp = T;

// Thermal postprocessing operations
        Tvisu = T - Tmelt;
        T_F = fvc::interpolate(T);
        gradT = fvc::snGrad(T);

        #if OPENFOAM >= 1906
        const MinMax<scalar> valTlimit = gMinMax(T);
        Info<< "Tmin= " << valTlimit.min() << nl
            << "Tmax= " << valTlimit.max() << endl;
        #else
        const scalar valTmin = gMin(T);
        const scalar valTmax = gMax(T);
        Info<< "Tmin= " << valTmin << nl
            << "Tmax= " << valTmax << endl;
        #endif


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
