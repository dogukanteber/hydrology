/*---------------------------------------------------------------------------*\
     Declaration of the variables used for the control of the picard Loops
\*---------------------------------------------------------------------------*/


const dictionary& picardDict = mesh.solutionDict().subDict("Picard");


//- multiplicating/dividing factor for time step changes
    const double tFact =
        picardDict.lookupOrDefault<double>("tFact", 1.);

//- maximum number of time step reductions for an iteration of the time loop
    const int nMaxCycle =
        picardDict.lookupOrDefault<int>("nMaxCycle", 1);


//- Picard loop parameters for Richards equation
//
//  - nIterPicard: maximum number of Picard loop iterations for the
//    linearisation of Richards equation
//  - precPicard: exit criterium of the Picard loop for the linearisation of
//    Richards equation
//  - stabilisationThresholdFlow: number of unchanged time steps with less than
//    three iterations of the Picard linearisation loop before of Richards
//    equation multiplying the time step by tFact (if the same conditions is
//    true for the heat transfer equation)
    const int nIterPicardFlow =
        picardDict.lookupOrDefault<int>("nIterPicardFlow", 1);
    const double precPicardFlow =
        picardDict.lookupOrDefault<double>("precPicardFlow", 1.);
    const int stabilisationThresholdFlow =
        picardDict.lookupOrDefault<int>("stabilisationThresholdFlow", 1);


//- Picard loop parameters for the heat transfer equation
//
//  - nIterPicardThermal: maximum number of Picard loop iterations for the
//    linearisation of the heat transfer equation
//  - precPicardThermal: exit criterium of the Picard loop for the
//    linearisation of the heat transfer equation
//  - stabilisationThresholdThermal: number of unchanged time steps with less
//    than three iterations of the Picard linearisation loop of the heat
//    transfer equation before multiplying the time step by tFact (if the same
//    conditions is true for Richards equation)
    const int nIterPicardThermal =
        picardDict.lookupOrDefault<int>("nIterPicardThermal", 1);
    const double precPicardThermal =
        picardDict.lookupOrDefault<double>("precPicardThermal", 1);
    const int stabilisationThresholdThermal =
        picardDict.lookupOrDefault<int>("stabilisationThresholdThermal", 1);





