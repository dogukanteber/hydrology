/*---------------------------------------------------------------------------*\
           Update some variables used to control of the Picard loops

Description
    Variables previously declared and initalized in declarePicardControls.H

Note
    This is a bit unwieldy, but we wish to react to runTime modifications
    of the parameters, but also remain 'sticky' between time loops
    if nothing has changed.
\*---------------------------------------------------------------------------*/

{
    constexpr int nIterPicardLimit = 3;

    bool emitMessage = false;

    if (nIterPicardFlow_old != nIterPicardFlow)
    {
        nIterPicardFlow_old = nIterPicardFlow;
        currentPicardFlow = nIterPicardFlow - nIterPicardLimit;
        stabilityCtrlFlow = 0;
        emitMessage = true;
    }

    if (emitMessage)
    {
        Info<< "Picard iterations. flow = " << nIterPicardFlow
            << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
