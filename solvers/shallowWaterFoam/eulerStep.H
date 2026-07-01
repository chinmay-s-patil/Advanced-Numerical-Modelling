// eulerStep.H -- Forward Euler for the LINEARIZED Shallow Water Equations
//
//   d(zeta)/dt  = -div(qFlux)
//   d(qFlux)/dt = f x perp(qFlux) - g*H*grad(zeta) - r*qFlux + windStress

// Wind gate: sum contributions from all active sources
const scalar t_ = runTime.value();
vector totalWindStress = vector::zero;
forAll(windSpeeds, i)
{
    if (t_ >= windT0s[i] && t_ <= windT1s[i])
    {
        totalWindStress += lambda.value() * windSpeeds[i] * windSpeeds[i] * windDirs[i];
    }
}
const dimensionedVector activeWind
(
    "activeWind",
    dimVelocity*dimVelocity,
    totalWindStress
);

// --- Continuity: d(zeta)/dt = -div(qFlux) ---
surfaceScalarField phiQ(fvc::interpolate(qFlux) & mesh.Sf());
volScalarField RHS_zeta_e = -fvc::div(phiQ);

// --- Coriolis ---
volVectorField cor
(
    IOobject("cor", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", qFlux.dimensions()/dimTime, vector::zero)
);
forAll(cor, ci)
{
    cor[ci] = f.value() * vector(qFlux[ci].y(), -qFlux[ci].x(), 0.0);
}
cor.correctBoundaryConditions();

// --- Linear friction: -r * qFlux (matches lecture equations) ---
volVectorField friction_e
(
    IOobject("friction_e", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    -r * qFlux
);

// --- Momentum RHS ---
volVectorField RHS_qFlux_e =
    cor
  - magg * H * fvc::grad(zeta)
  + friction_e
  + activeWind;

// --- Euler update ---
zeta  += runTime.deltaT() * RHS_zeta_e;
qFlux += runTime.deltaT() * RHS_qFlux_e;

// --- Mask land cells (H = 0) ---
forAll(H, cellI)
{
    if (H[cellI] < SMALL)
    {
        qFlux[cellI] = vector::zero;
        zeta[cellI]  = 0.0;
    }
}

zeta.correctBoundaryConditions();
qFlux.correctBoundaryConditions();
