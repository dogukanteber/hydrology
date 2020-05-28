/*---------------------------------------------------------------------------*\
                  Handling of adaptive time steping

Description
    Reset the timestep to maintain good convergence behaviours in non
      linearities handling by picard Loops.
    Reduction of time-step is immediate, but increase is damped to avoid
      unstable oscillations.

\*---------------------------------------------------------------------------*/

scalar deltaTFact = 1.;


if (adjustTimeStep)
{
    deltaTFact = 1.;


    if
    (
        (currentPicardFlow > nIterPicardFlow -2) ||
        (currentPicardThermal > nIterPicardThermal -2)
    )
    {
        deltaTFact = 1./tFact;
    }


    if (currentPicardFlow >= 3)
    {
        scf = 0;
    }


    if (currentPicardThermal >= 3)
    {
        scth = 0;
    }


    if (currentPicardFlow < 3)
    {
        scf = scf+1;
    }


    if (currentPicardThermal < 3)
    {
        scth = scth+1;
    }

    if
    (
        (
            (currentPicardFlow < 3) &&
            (scf >= stabilisationThresholdFlow)
        )
        &&
        (
            (currentPicardThermal < 3) &&
            (scth >= stabilisationThresholdThermal)
        )
    )
    {
        deltaTFact = tFact;
        scf = 0;
        scth = 0;
    }


    if (deltaTFact != 1.)
    {
        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        );
    }


    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}


// ************************************************************************* //
