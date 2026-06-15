Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/constructing bathymetry field zb\n" << endl;
volScalarField zb
(
    IOobject
    (
        "zb",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zb", dimLength, 0.0)
);

Info << "Reading SWEDict\n" << endl;
IOdictionary SWEDict
(
    IOobject
    (
        "SWEDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const word type
(
    SWEDict.get<word>("type")
);

Info<< "Reading g\n" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar magg("magg", mag(g));

//  Coriolis parameter 
// Default: f = 2 * 7.2921e-5 * sin(60 deg) = 1.2585e-4 s^-1  (60 N)
const dimensionedScalar f
(
    "f",
    dimless/dimTime,
    SWEDict.lookupOrDefault<scalar>("f", 1.2585e-4)
);

//  Bottom friction coefficient 
// Drag term: -r * hU * |hU| / h^2
const dimensionedScalar r
(
    "r",
    dimless,
    SWEDict.lookupOrDefault<scalar>("r", 0.003)
);

//  Wind stress factor 
const dimensionedScalar lambda
(
    "lambda",
    dimless,
    SWEDict.lookupOrDefault<scalar>("lambda", 3.2e-6)
);

//  Source (wind) subdict 
//
// In SWEDict, add:
//
//   source
//   {
//       speed       10;          // wind speed [m/s]
//       direction   (1 0 0);     // unit vector: East =(1,0,0)
//                                //              North=(0,1,0)
//                                //              West =(-1,0,0)
//                                //              South=(0,-1,0)
//       duration    (0 1e9);     // (tStart tEnd) [s] — wind active window
//   }

const dictionary& sourceDict = SWEDict.subDict("source");

const scalar windSpeed = sourceDict.get<scalar>("speed");

const vector windDirRaw = sourceDict.get<vector>("direction");
const vector windDirUnit = windDirRaw / (mag(windDirRaw) + SMALL);

// Duration window
const scalarList windDurationList =
    sourceDict.get<scalarList>("duration");

const scalar windT0 = windDurationList[0];   // start time [s]
const scalar windT1 = windDurationList[1];   // end   time [s]

// Pre-compute constant wind stress vector: lambda * W * |W|
// Units: [m^2/s^2] — added directly to the hU momentum equation RHS
const dimensionedVector windStress
(
    "windStress",
    dimVelocity * dimVelocity,
    lambda.value() * windSpeed * windSpeed * windDirUnit
);

Info<< "Wind: speed=" << windSpeed << " m/s"
    << "  direction=" << windDirUnit
    << "  active=[" << windT0 << ", " << windT1 << "] s" << endl;

//  Conserved momentum field hU = h*U 
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h*U
);
