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

//- Temperature field [K] at cell centers
//  primary variable of resolution for thermal resolution
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
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

//- Temperature field [K] - face centers interpolated
    Info<< "Reading field T_F\n" << endl;
    surfaceScalarField T_F
    (
        IOobject
        (
            "T_F",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Temperature gradient field [K.m-1] - projected on faces
    Info<< "Reading field gradT\n" << endl;
    surfaceScalarField gradT
    (
        IOobject
        (
            "gradT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Linearization error associated with thermal resolution [K]
    Info<< "Reading field errT\n" << endl;
    volScalarField errT
    (
        IOobject
        (
            "errT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Temperature field in celsius degree [C]
    Info<< "Reading field Tvisu\n" << endl;
    volScalarField Tvisu
    (
        IOobject
        (
            "Tvisu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Ice volumetric content [-]
    Info<< "Reading field thetag\n" << endl;
    volScalarField thetag
    (
        IOobject
        (
            "thetag",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//- Liquid water content [-]
    Info<< "Reading field thetal\n" << endl;
    volScalarField thetal
    (
        IOobject
        (
            "thetal",
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh
    );

//- Soil thermal conductivity field [kg.m.s-3.K-1]
    Info<< "Reading field Kthsoil\n" << endl;
    volScalarField Kthsoil
    (
        IOobject
        (
            "Kthsoil",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

//- Soil heat capacity field [kg.m-1.s-2.K-1]
    Info<< "Reading field Cthsoil\n" << endl;
    volScalarField Cthsoil
    (
        IOobject
        (
            "Cthsoil",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

//- Temperature of start of the freezing [K]
    Info<< "Reading field Tmelt\n" << endl;
    volScalarField Tmelt
    (
        IOobject
        (
            "Tmelt",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

//- Minimal frozen relative hydraulic conductivity [-]
    Info<< "Reading field KrThMin\n" << endl;
    volScalarField KrThMin
    (
        IOobject
        (
            "KrThMin",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

//- Freezing impedance factor for frozen relative hydraulic conductivity [-]
    Info<< "Reading field omega\n" << endl;
    volScalarField omega
    (
        IOobject
        (
            "omega",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

//- Scaling parameter of the freezing curve of the soil [K]
    Info<< "Reading field W\n" << endl;
    volScalarField W
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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

//- Liquid water thermal conductivity [kg.m.s-3.K-1]
    dimensionedScalar Kthw
    (
        transportProperties.lookup("Kthw")
    );

//- Air thermal conductivity [kg.m.s-3.K-1]
    dimensionedScalar Kthair
    (
        transportProperties.lookup("Kthair")
    );

//- Ice thermal conductivity [kg.m.s-3.K-1]
    dimensionedScalar Kthice
    (
        transportProperties.lookup("Kthice")
    );

//- Liquid water heat capacity [kg.m-1.s-2.K-1]
    dimensionedScalar Cthw
    (
        transportProperties.lookup("Cthw")
    );

//- Air heat capacity [kg.m-1.s-2.K-1]
    dimensionedScalar Cthair
    (
        transportProperties.lookup("Cthair")
    );

//- Ice heat capacity [kg.m-1.s-2.K-1]
   dimensionedScalar Cthice
    (
        transportProperties.lookup("Cthice")
    );

//- Latent of fusion of ice [kg.m-1.s-2]
   dimensionedScalar L
    (
        transportProperties.lookup("L")
    );

//- Dimensionnning factor for the computation of apparent thermal conductivity
//  allways constant equal to 1 [kg.m.s-3.K-1]
     dimensionedScalar Kthdim
     (
         transportProperties.lookup("Kthdim")
     );

//- Dimensionnning factor for the computation of actual evapotranspiration AET
//  allways constant equal to 1 [s]
    dimensionedScalar tdim
    (
        transportProperties.lookup("tdim")
    );

//- Dimensionnning factor for the phase changes tests
//  allways constant equal to 1 [K]
    dimensionedScalar T1
    (
        transportProperties.lookup("T1")
    );

//- Maximum temperature (forced bound) [K]
    dimensionedScalar Tmaxc
    (
        transportProperties.lookup("Tmaxc")
    );

//- Minimum temperature (forced bound) [K]
    dimensionedScalar Tminc
    (
        transportProperties.lookup("Tminc")
    );
