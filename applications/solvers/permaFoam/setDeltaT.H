/*---------------------------------------------------------------------------*\
                  Handling of adaptive time steping

Description
    Reset the timestep to maintain good convergence behaviours in
    non-linearities handling by Picard Loops.

    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    constexpr int nIterPicardLimit = 3;

    scalar deltaTFactor = 1;

    if
    (
        (currentPicardFlow > nIterPicardFlow - nIterPicardLimit + 1)
     || (currentPicardThermal > nIterPicardThermal - nIterPicardLimit + 1)
    )
    {
        deltaTFactor = 1./tFact;
    }


    if (currentPicardFlow < nIterPicardLimit)
    {
        ++stabilityCtrlFlow;
    }
    else
    {
        stabilityCtrlFlow = 0;  // Reset
    }


    if (currentPicardThermal < nIterPicardLimit)
    {
        ++stabilityCtrlThermal;
    }
    else
    {
        stabilityCtrlThermal = 0;  // Reset
    }

    if
    (
        (
            (currentPicardFlow < nIterPicardLimit)
         && (stabilityCtrlFlow >= stabilityThresholdFlow)
        )
     &&
        (
            (currentPicardThermal < nIterPicardLimit)
         && (stabilityCtrlThermal >= stabilityThresholdThermal)
        )
    )
    {
        deltaTFactor = tFact;
        // Reset stability controls
        stabilityCtrlFlow = 0;
        stabilityCtrlThermal = 0;
    }


    if (!equal(deltaTFactor, 1))
    {
        runTime.setDeltaT
        (
            min
            (
                deltaTFactor*runTime.deltaTValue(),
                maxDeltaT
            )
        );
    }
}

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;


// ************************************************************************* //
