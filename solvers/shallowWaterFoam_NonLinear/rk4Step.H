// rk4Step.H
// Classical RK4 with Coriolis, bottom friction, and wind stress.
//
// Source term evaluated at each stage:
//   S(h, hU) = + f*perp(hU)              Coriolis
//              - r * hU * |hU| / h^2     friction
//              + activeWind              wind stress (time-gated)

// Wind on/off at current time (same for all 4 stages — time-centred)
const scalar t_ = runTime.value();
const dimensionedVector activeWind
(
    "activeWind",
    dimVelocity*dimVelocity,
    (t_ >= windT0 && t_ <= windT1) ? windStress.value() : vector::zero
);

// Dry-cell threshold: must match the initial "dry" depth in 0/h
const dimensionedScalar hDry("hDry", dimLength, 0.05);

// Maximum allowed velocity magnitude [m/s] — physical bound for lake flows
const scalar Umax = 50.0;

//  Helper: clamp h, zero momentum in dry cells, limit velocity magnitude
auto clampFields = [&]()
{
    // Floor depth
    h = max(h, hDry);
    h.correctBoundaryConditions();

    // Recover velocity
    U = hU / h;

    // Zero velocity in dry cells (h near hDry)
    volScalarField wetFrac
    (
        Foam::max
        (
            (h - hDry) / dimensionedScalar("dh", dimLength, 1.0),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    wetFrac = Foam::min(wetFrac, dimensionedScalar("one", dimless, 1.0));
    U *= wetFrac;

    // Clamp velocity magnitude
    volScalarField magU(mag(U));
    volScalarField scale
    (
        Foam::min
        (
            magU / dimensionedScalar("Umax", dimVelocity, Umax) + dimensionedScalar("small", dimless, SMALL),
            dimensionedScalar("one", dimless, 1.0)
        )
    );
    // scale = min(|U|/Umax, 1). If |U| > Umax, we want to reduce by Umax/|U|
    // So: U_clamped = U * min(1, Umax / |U|)
    volScalarField clampFactor
    (
        Foam::min
        (
            dimensionedScalar("Umax", dimVelocity, Umax) / (magU + dimensionedScalar("small", dimVelocity, SMALL)),
            dimensionedScalar("one", dimless, 1.0)
        )
    );
    U *= clampFactor;

    U.correctBoundaryConditions();
    hU = h * U;
};

//  Helper: compute hU source given a (h_, hU_) pair
auto sweSourceHU = [&](
    const volScalarField& h_,
    const volVectorField& hU_) -> volVectorField
{
    // Coriolis
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

    // Friction
    volScalarField fMag
    (
        mag(hU_) /
        (h_*h_ + dimensionedScalar("small", dimArea, SMALL))
    );

    volScalarField fMagLim
    (
        min
        (
            fMag,
            1.0 / (r*runTime.deltaT() + dimensionedScalar("smallTime", dimTime, SMALL))
        )
    );

    // Wind only on wet cells
    volScalarField windMask
    (
        pos0(h_ - 2.0 * hDry)
    );

    return cor - r*fMagLim*hU_ + windMask*activeWind;
};

//  Snapshot
h0  = h;
hU0 = hU;

//  Stage 1
{
    surfaceScalarField phi1(fvc::interpolate(hU) & mesh.Sf());

    k1_h  = -fvc::div(phi1);
    k1_hU = -fvc::div(phi1, U)
           - (magg/2)*fvc::grad(h*h)
           - magg*h*fvc::grad(zb)
           + sweSourceHU(h, hU)
           + h * nuArt * fvc::laplacian(U);
}

//  Stage 2
h  = h0  + (runTime.deltaT()/2)*k1_h;
hU = hU0 + (runTime.deltaT()/2)*k1_hU;
clampFields();
{
    surfaceScalarField phi2(fvc::interpolate(hU) & mesh.Sf());

    k2_h  = -fvc::div(phi2);
    k2_hU = -fvc::div(phi2, U)
           - (magg/2)*fvc::grad(h*h)
           - magg*h*fvc::grad(zb)
           + sweSourceHU(h, hU)
           + h * nuArt * fvc::laplacian(U);
}

//  Stage 3
h  = h0  + (runTime.deltaT()/2)*k2_h;
hU = hU0 + (runTime.deltaT()/2)*k2_hU;
clampFields();
{
    surfaceScalarField phi3(fvc::interpolate(hU) & mesh.Sf());

    k3_h  = -fvc::div(phi3);
    k3_hU = -fvc::div(phi3, U)
           - (magg/2)*fvc::grad(h*h)
           - magg*h*fvc::grad(zb)
           + sweSourceHU(h, hU)
           + h * nuArt * fvc::laplacian(U);
}

//  Stage 4
h  = h0  + (runTime.deltaT()/2)*k3_h;
hU = hU0 + (runTime.deltaT()/2)*k3_hU;
clampFields();
{
    surfaceScalarField phi4(fvc::interpolate(hU) & mesh.Sf());

    k4_h  = -fvc::div(phi4);
    k4_hU = -fvc::div(phi4, U)
           - (magg/2)*fvc::grad(h*h)
           - magg*h*fvc::grad(zb)
           + sweSourceHU(h, hU)
           + h * nuArt * fvc::laplacian(U);
}

//  Final update
h  = h0  + (runTime.deltaT()/6)*(k1_h  + 2*k2_h  + 2*k3_h  + k4_h);
hU = hU0 + (runTime.deltaT()/6)*(k1_hU + 2*k2_hU + 2*k3_hU + k4_hU);

clampFields();
