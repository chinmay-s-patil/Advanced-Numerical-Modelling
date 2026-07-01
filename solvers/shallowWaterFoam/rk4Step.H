// rk4Step.H -- Classical RK4 for the LINEARIZED Shallow Water Equations
//
// Linearized SWE (lecture equations):
//   d(zeta)/dt  = -div(qFlux)
//   d(qFlux)/dt = f x perp(qFlux) - g*H*grad(zeta) - r*qFlux + windStress
//
// No advection, no nonlinear terms. Simple and clean.

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

// --- RHS evaluation lambda ---
// Returns (d(zeta)/dt, d(qFlux)/dt) for a given state (zeta_, qFlux_)

auto RHS_zeta = [&](
    const volScalarField& /*zeta_*/,
    const volVectorField& qFlux_) -> volScalarField
{
    // Continuity: d(zeta)/dt = -div(qFlux)
    surfaceScalarField phiQ(fvc::interpolate(qFlux_) & mesh.Sf());
    return -fvc::div(phiQ);
};

auto RHS_qFlux = [&](
    const volScalarField& zeta_,
    const volVectorField& qFlux_) -> volVectorField
{
    // Coriolis: f x perp(qFlux) = f * (qFlux.y, -qFlux.x, 0)
    volVectorField cor
    (
        IOobject("cor", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("zero", qFlux_.dimensions()/dimTime, vector::zero)
    );
    forAll(cor, ci)
    {
        cor[ci] = f.value() * vector(qFlux_[ci].y(), -qFlux_[ci].x(), 0.0);
    }
    cor.correctBoundaryConditions();

    // Pressure gradient: -g*H*grad(zeta)
    volVectorField pressGrad = -magg * H * fvc::grad(zeta_);

    // Linear friction: -r * qFlux (matches lecture equations)
    volVectorField friction
    (
        IOobject("friction", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        -r * qFlux_
    );

    return cor + pressGrad + friction + activeWind;
};

// --- Snapshot current state ---
zeta0  = zeta;
qFlux0 = qFlux;

// --- Stage 1 ---
k1_zeta  = RHS_zeta(zeta, qFlux);
k1_qFlux = RHS_qFlux(zeta, qFlux);

// --- Stage 2 ---
zeta  = zeta0  + (runTime.deltaT()/2)*k1_zeta;
qFlux = qFlux0 + (runTime.deltaT()/2)*k1_qFlux;
zeta.correctBoundaryConditions();
qFlux.correctBoundaryConditions();

k2_zeta  = RHS_zeta(zeta, qFlux);
k2_qFlux = RHS_qFlux(zeta, qFlux);

// --- Stage 3 ---
zeta  = zeta0  + (runTime.deltaT()/2)*k2_zeta;
qFlux = qFlux0 + (runTime.deltaT()/2)*k2_qFlux;
zeta.correctBoundaryConditions();
qFlux.correctBoundaryConditions();

k3_zeta  = RHS_zeta(zeta, qFlux);
k3_qFlux = RHS_qFlux(zeta, qFlux);

// --- Stage 4 ---
zeta  = zeta0  + runTime.deltaT()*k3_zeta;
qFlux = qFlux0 + runTime.deltaT()*k3_qFlux;
zeta.correctBoundaryConditions();
qFlux.correctBoundaryConditions();

k4_zeta  = RHS_zeta(zeta, qFlux);
k4_qFlux = RHS_qFlux(zeta, qFlux);

// --- Final update ---
zeta  = zeta0  + (runTime.deltaT()/6)*(k1_zeta  + 2*k2_zeta  + 2*k3_zeta  + k4_zeta);
qFlux = qFlux0 + (runTime.deltaT()/6)*(k1_qFlux + 2*k2_qFlux + 2*k3_qFlux + k4_qFlux);

// --- Mask land cells (H = 0) ---
// Land cells must have zero transport and zero free-surface perturbation
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
