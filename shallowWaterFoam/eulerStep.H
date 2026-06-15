// eulerStep.H

// Setting phi from hU (dot) Surface Area Vector
surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(hU) & mesh.Sf()
);

// Solving the h Equation (Continuity)
solve
(
    fvm::ddt(h) + fvc::div(phi)
);

// Solving Momentum Equation
solve
(
    fvm::ddt(hU)
  + fvc::div(phi, U)
  + (magg/2)*fvc::grad(h*h)
  + magg*h*fvc::grad(zb)
);

// Recover U
U = hU/h;
