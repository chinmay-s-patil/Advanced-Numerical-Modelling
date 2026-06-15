// cnStep.H
// Explicit Crank-Nicolson (predictor-corrector) with Coriolis, friction, wind.
//
// Same source term pattern as eulerStep / rk4Step:
//   S(h, hU) = f*perp(hU) - r*|hU|*hU/h^2 + activeWind

// Wind gate
const scalar t_ = runTime.value();
const dimensionedVector activeWind
(
    "activeWind",
    dimVelocity*dimVelocity,
    (t_ >= windT0 && t_ <= windT1) ? windStress.value() : vector::zero
);

// Reusable source-term lambda (same pattern as rk4Step)
auto sweSourceHU = [&](
    const volScalarField& h_,
    const volVectorField& hU_) -> volVectorField
{
    volVectorField cor
    (
        IOobject("cor", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("zero", hU_.dimensions()/dimTime, vector::zero)
    );
    forAll(cor, ci)
    {
        cor[ci] = f.value() * vector(hU_[ci].y(), -hU_[ci].x(), 0.0);
    }
    cor.correctBoundaryConditions();

    volScalarField fMag
    (
        mag(hU_) /
        (h_*h_ + dimensionedScalar("small", dimArea, SMALL))
    );

    return cor - r*fMag*hU_ + activeWind;
};

// ── Step 1: RHS at current state (n) ─────────────────────────────────────
surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(hU) & mesh.Sf()
);

RHS_h  = fvc::div(-phi);
RHS_hU = -fvc::div(phi, U)
        - (magg/2)*fvc::grad(h*h)
        - magg*h*fvc::grad(zb)
        + sweSourceHU(h, hU);

// ── Step 2: Euler predictor ───────────────────────────────────────────────
hStar  = h  + runTime.deltaT()*RHS_h;
hUStar = hU + runTime.deltaT()*RHS_hU;
UStar  = hUStar/hStar;

// ── Step 3: RHS at star state ─────────────────────────────────────────────
surfaceScalarField phiStar
(
    IOobject("phiStar", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(hUStar) & mesh.Sf()
);

RHS_h_Star  = fvc::div(-phiStar);
RHS_hU_Star = -fvc::div(phiStar, UStar)
             - (magg/2)*fvc::grad(hStar*hStar)
             - magg*hStar*fvc::grad(zb)
             + sweSourceHU(hStar, hUStar);

// ── Step 4: Corrector average ─────────────────────────────────────────────
h  = h  + (runTime.deltaT()/2)*(RHS_h  + RHS_h_Star);
hU = hU + (runTime.deltaT()/2)*(RHS_hU + RHS_hU_Star);

// ── Step 5: Recover U ─────────────────────────────────────────────────────
U = hU/h;
