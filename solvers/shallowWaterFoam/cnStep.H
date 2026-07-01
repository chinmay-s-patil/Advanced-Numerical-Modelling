// cnStep.H -- Explicit Crank-Nicolson (predictor-corrector) for LINEARIZED SWE
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

// --- RHS evaluation helper ---
auto RHS_zeta_cn = [&](
    const volScalarField& /*zeta_*/,
    const volVectorField& qFlux_) -> volScalarField
{
    surfaceScalarField phiQ(fvc::interpolate(qFlux_) & mesh.Sf());
    return -fvc::div(phiQ);
};

auto RHS_qFlux_cn = [&](
    const volScalarField& zeta_,
    const volVectorField& qFlux_) -> volVectorField
{
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

    // Linear friction: -r * qFlux (matches lecture equations)
    volVectorField friction_cn
    (
        IOobject("friction_cn", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        -r * qFlux_
    );

    return cor - magg * H * fvc::grad(zeta_) + friction_cn + activeWind;
};

// Step 1: RHS at current state (n)
volScalarField  Rz_n  = RHS_zeta_cn(zeta, qFlux);
volVectorField  Rq_n  = RHS_qFlux_cn(zeta, qFlux);

// Step 2: Euler predictor
volScalarField  zetaStar  = zeta  + runTime.deltaT()*Rz_n;
volVectorField  qFluxStar = qFlux + runTime.deltaT()*Rq_n;
zetaStar.correctBoundaryConditions();
qFluxStar.correctBoundaryConditions();

// Step 3: RHS at star state
volScalarField  Rz_star = RHS_zeta_cn(zetaStar, qFluxStar);
volVectorField  Rq_star = RHS_qFlux_cn(zetaStar, qFluxStar);

// Step 4: Corrector average
zeta  += (runTime.deltaT()/2)*(Rz_n  + Rz_star);
qFlux += (runTime.deltaT()/2)*(Rq_n  + Rq_star);

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
