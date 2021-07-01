/*---------------------------------------------------------------------------*\
                    Declare some variables for Picard loops
\*---------------------------------------------------------------------------*/

// Variables for adaptative time stepping procedure

int nIterPicardFlow_old = -1;
int nIterPicardThermal_old = -1;

int currentPicardFlow = 0;
int currentPicardThermal = 0;

// Stabilisation controls
int stabilityCtrlFlow = 0;
int stabilityCtrlThermal = 0;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
