// createFields.H -- Linearized Shallow Water Equations
// Prognostic variables:
//   zeta  -- free-surface height perturbation [m]
//   qFlux -- depth-integrated transport vector (U, V, 0) [m^2/s]
// Static field:
//   H     -- bathymetric depth [m]  (read from file, never changes)

Info<< "Reading bathymetric depth field H\n" << endl;
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

Info<< "Reading/initialising free-surface perturbation zeta\n" << endl;
volScalarField zeta
(
    IOobject
    (
        "zeta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/initialising depth-integrated transport qFlux\n" << endl;
volVectorField qFlux
(
    IOobject
    (
        "qFlux",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

// ---------------------------------------------------------------------------
// Physical parameters
// ---------------------------------------------------------------------------

// Coriolis parameter [s^-1]
// f = 2 * Omega * sin(latitude)
// At 60 N:  f = 2 * 7.2921e-5 * sin(60) = 1.2585e-4
const dimensionedScalar f
(
    "f",
    dimless/dimTime,
    SWEDict.lookupOrDefault<scalar>("f", 1.2585e-4)
);

// Bottom friction coefficient (linear drag)
// Linearized friction: -r * qFlux
const dimensionedScalar r
(
    "r",
    dimless/dimTime,
    SWEDict.lookupOrDefault<scalar>("r", 0.003)
);

// Wind stress factor lambda
const dimensionedScalar lambda
(
    "lambda",
    dimless,
    SWEDict.lookupOrDefault<scalar>("lambda", 3.2e-6)
);

// ---------------------------------------------------------------------------
// Source (wind) subdict -- supports multiple simultaneous sources
// ---------------------------------------------------------------------------
//
// In SWEDict, use either single source (backward compatible):
//
//   source
//   {
//       speed       10;
//       direction   (1 0 0);
//       duration    (0 1e9);
//   }
//
// Or multiple sources:
//
//   source
//   {
//       easterly
//       {
//           speed       10;
//           direction   (1 0 0);
//           duration    (0 50);
//       }
//       northerly
//       {
//           speed       5;
//           direction   (0 1 0);
//           duration    (0 50);
//       }
//   }

const dictionary& sourceDict = SWEDict.subDict("source");

List<scalar> windSpeeds(0);
List<vector> windDirs(0);
List<scalar> windT0s(0);
List<scalar> windT1s(0);

// Check if sourceDict has sub-dictionaries (multi-source) or direct entries (single-source)
bool hasSubDicts = false;
forAllConstIter(dictionary, sourceDict, iter)
{
    if (iter().isDict())
    {
        hasSubDicts = true;
        break;
    }
}

if (hasSubDicts)
{
    // Multi-source mode
    forAllConstIter(dictionary, sourceDict, iter)
    {
        const entry& e = iter();
        if (e.isDict())
        {
            const dictionary& subDict = e.dict();

            scalar speed = subDict.get<scalar>("speed");
            vector dirRaw = subDict.get<vector>("direction");
            vector dirUnit = dirRaw / (mag(dirRaw) + SMALL);
            scalarList dur = subDict.get<scalarList>("duration");

            label n = windSpeeds.size();
            windSpeeds.setSize(n+1);
            windDirs.setSize(n+1);
            windT0s.setSize(n+1);
            windT1s.setSize(n+1);

            windSpeeds[n] = speed;
            windDirs[n] = dirUnit;
            windT0s[n] = dur[0];
            windT1s[n] = dur[1];

            Info<< "Wind source '" << e.keyword() << "': speed=" << speed 
                << " m/s, direction=" << dirUnit 
                << ", active=[" << dur[0] << ", " << dur[1] << "] s" << endl;
        }
    }
}
else
{
    // Single-source mode (backward compatible)
    scalar speed = sourceDict.get<scalar>("speed");
    vector dirRaw = sourceDict.get<vector>("direction");
    vector dirUnit = dirRaw / (mag(dirRaw) + SMALL);
    scalarList dur = sourceDict.get<scalarList>("duration");

    windSpeeds.setSize(1);
    windDirs.setSize(1);
    windT0s.setSize(1);
    windT1s.setSize(1);

    windSpeeds[0] = speed;
    windDirs[0] = dirUnit;
    windT0s[0] = dur[0];
    windT1s[0] = dur[1];

    Info<< "Wind source (single): speed=" << speed 
        << " m/s, direction=" << dirUnit 
        << ", active=[" << dur[0] << ", " << dur[1] << "] s" << endl;
}

Info<< "Total wind sources: " << windSpeeds.size() << endl;
