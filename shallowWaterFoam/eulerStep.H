// eulerStep.H
// Explicit Forward Euler with Coriolis, bottom friction, and wind stress.
//
// Momentum RHS extra terms (from lecture eq. 1-2):
//   Coriolis    : + f * perp(hU)          perp(hU) = (hU.y, -hU.x, 0)
//   Friction    : - r * hU * |hU| / h^2
//   Wind stress : + windStress             (active only within duration window)

// Wind on/off switch for current time
const scalar t_ = runTime.value();
const dimensionedVector activeWind
(
    "activeWind",
    dimVelocity*dimVelocity,
    (t_ >= windT0 && t_ <= windT1) ? windStress.value() : vector::zero
);

// Coriolis: rotate hU by 90 deg in the horizontal plane
// perp(hU) = (hU.y, -hU.x, 0)
volVectorField coriolisForce
(
    "coriolisForce",
    f * volVectorField
    (
        IOobject("tmp", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("zero", hU.dimensions()/dimTime, vector::zero)
    )
);
forAll(coriolisForce, ci)
{
    coriolisForce[ci] = f.value() *
        vector(hU[ci].y(), -hU[ci].x(), 0.0);
}
coriolisForce.correctBoundaryConditions();

// Bottom friction magnitude: |hU| / h^2
volScalarField frictionMag
(
    "frictionMag",
    mag(hU) / (h*h + dimensionedScalar("small", dimArea, SMALL))
);

// Flux
surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(hU) & mesh.Sf()
);

// Continuity
solve
(
    fvm::ddt(h) + fvc::div(phi)
);

// Momentum
solve
(
    fvm::ddt(hU)
  + fvc::div(phi, U)
  + (magg/2)*fvc::grad(h*h)
  + magg*h*fvc::grad(zb)
  - coriolisForce                          // Coriolis
  + r*frictionMag*hU                      // bottom friction
  - activeWind                            // wind stress
);

// Recover U
U = hU/h;
