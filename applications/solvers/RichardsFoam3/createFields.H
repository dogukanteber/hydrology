/*---------------------------------------------------------------------------*\
                Declarations of the fields of primary variables
\*---------------------------------------------------------------------------*/


//- Total water pressure (pressure head/matrix potential) [m] at cell centers
//  primary variable of resolution for Richards equation
    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Piezometric head [m] at cell centers
    Info<< "Reading field H\n" << endl;
    volScalarField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


/*---------------------------------------------------------------------------*\
                 Declarations of the fields of derived variables
\*---------------------------------------------------------------------------*/


//- Water pressure field [m] - face centers interpolated
    Info<< "Reading field psi_F\n" << endl;
    surfaceScalarField psi_F
    (
        IOobject
        (
            "psi_F",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Water velocity field normal to faces [m/s]
    Info<< "Reading field U\n" << endl;
    surfaceScalarField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Velocity vector field at the center of cells [m/s]
    Info<< "Reading field Uvol\n" << endl;
    volVectorField Uvol
    (
        IOobject
        (
            "Uvol",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Total water volumetric content [-]
    Info<< "Reading field theta\n" << endl;
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Actual evapotranspiration [s-1]
    Info<< "Reading field AET\n" << endl;
    volScalarField AET
    (
        IOobject
        (
            "AET",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Field of linearization error for the resolution of water transfers [m]
    Info<< "Reading field err\n" << endl;
    volScalarField err
    (
        IOobject
        (
            "err",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

/*---------------------------------------------------------------------------*\
     Declarations of the fields of input data (medium properties, forcings)
\*---------------------------------------------------------------------------*/


//- Saturated hydraulic conductivity field [m/s]
    Info<< "Reading field K\n" << endl;
    volScalarField K
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Inverse of the capillary length [m-1]
    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- van Genuchten exponent [-]
    Info<< "Reading field n\n" << endl;
    volScalarField n
    (
        IOobject
        (
            "n",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Saturation water content [-]
    Info<< "Reading field thetas\n" << endl;
    volScalarField thetas
    (
        IOobject
        (
            "thetas",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Residual water content [-]
    Info<< "Reading field thetar\n" << endl;
    volScalarField thetar
    (
        IOobject
        (
            "thetar",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Water storage coefficient [m-1]
    Info<< "Reading field S\n" << endl;
    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Potential evapotranspiration [s-1]
//  in each cell of the uptake zone:
//  PET = PET flux [m.s-1] * total top surface [m2] / uptake zone volume [m3]
//  anywhere else (i.e., deeper than the uptake zone): PET = 0
    Info<< "Reading field PET\n" << endl;
    volScalarField PET
    (
        IOobject
        (
            "PET",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Water content at wilting point [-]
    Info<< "Reading field thetaWP\n" << endl;
    volScalarField thetaWP
    (
        IOobject
        (
            "thetaWP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

/*---------------------------------------------------------------------------*\
     Geometrical fields
\*---------------------------------------------------------------------------*/

//- Vertical ascendant unit vector [-]
    Info<< "Reading field vuz\n" << endl;
    volVectorField vuz
    (
        IOobject
        (
            "vuz",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Unit metric scalar field [m]
    Info<< "Reading field usf\n" << endl;
    volScalarField usf
    (
        IOobject
        (
            "usf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


/*---------------------------------------------------------------------------*\
                        Reading of scalar input data
\*---------------------------------------------------------------------------*/


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


//- Dimensioning factor for the computation of actual evapotranspiration AET
//  always constant equal to 1 [s]
    const dimensionedScalar tdim("tdim", transportProperties);

// ************************************************************************* //
