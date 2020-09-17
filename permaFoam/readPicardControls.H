/*---------------------------------------------------------------------------*\
                 Variables used to control of the Picard loops
\*---------------------------------------------------------------------------*/
// Note
//     Uses dictionary lookupOrDefault() instead of getOrDefault()
//     to ensure compilation with OpenFOAM-v1812 and older.

const dictionary& picardDict = mesh.solutionDict().subDict("Picard");


//- Multiplication/division factor for time step changes
const scalar tFact =
    picardDict.lookupOrDefault<scalar>("tFact", 1);

//- Maximum number of time step reductions for an iteration of the time loop
const int nMaxCycle =
    picardDict.lookupOrDefault<int>("nMaxCycle", 1);


//- Picard loop parameters for Richards equation
//
//  - nIterPicardFlow:
//    maximum number of Picard loop iterations for the linearisation of
//    Richards equation
//  - precPicardFlow:
//    exit criterium of the Picard loop for the linearisation of
//    Richards equation
//  - stabilisationThresholdFlow:
//    number of unchanged time steps with less than three iterations of
//    the Picard linearisation loop
//    Richards equation before multiplying the time step by tFact
//    (if the same conditions is true for the heat transfer equation)
const int nIterPicardFlow =
    picardDict.lookupOrDefault<int>("nIterPicardFlow", 10);

const scalar precPicardFlow =
    picardDict.lookupOrDefault<scalar>("precPicardFlow", 1e-3);

const int stabilisationThresholdFlow =
    picardDict.lookupOrDefault<int>("stabilisationThresholdFlow", 50);


//- Picard loop parameters for the heat transfer equation
//
//  - nIterPicardThermal:
//    maximum number of Picard loop iterations for the linearisation of
//    the heat transfer equation
//  - precPicardThermal:
//    exit criterium of the Picard loop for the linearisation of
//    the heat transfer equation
//  - stabilisationThresholdThermal:
//    number of unchanged time steps with less than three iterations of
//    the Picard linearisation loop of
//    the heat transfer equation before multiplying the time step by tFact
//    (if the same conditions is true for Richards equation)
//
//    Use flow equation values for defaults
const int nIterPicardThermal =
    picardDict.lookupOrDefault<int>
    (
        "nIterPicardThermal",
        nIterPicardFlow
    );

const scalar precPicardThermal =
    picardDict.lookupOrDefault<scalar>
    (
        "precPicardThermal",
        precPicardFlow
    );

const int stabilisationThresholdThermal =
    picardDict.lookupOrDefault<int>
    (
        "stabilisationThresholdThermal",
        stabilisationThresholdFlow
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
