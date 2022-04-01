Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE 
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ), 
    mesh
);

Info<< "Reading field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

// ~~~~~~~~~~ transportProperties:
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// ~~~~~~~~~~ nu (read from):
dimensionedScalar nu
(   "nu",
    transportProperties.lookup("nu")
);

// ~~~~~~~~~~ etta (read from):
dimensionedScalar etta
    (
    "etta", 
    transportProperties.lookup("etta")
    );

// ~~~~~~~~~~ gamma:
Info<< "Reading field gamma\n" << endl;
volScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// ~~~~~~~~~~ DFx/Dgma
volScalarField DFx_Dgma
(
    IOobject
    (
        "DFx_Dgma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DFx_Dgma",  dimensionSet(0,0,0,0,0,0,0), 0.0)
);

// ~~~~~~~~~~ DFy/Dgma
volScalarField DFy_Dgma
(
    IOobject
    (
        "DFy_Dgma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DFy_Dgma",  dimensionSet(0,0,0,0,0,0,0), 0.0)
);

// ~~~~~~~~~~ p2 field for pEqn.solve():
volScalarField p2
(
    IOobject
    (
        "p2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p2",  dimensionSet(0,2,-2,0,0,0,0), scalar(0.0))
);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END

