volScalarField h0("h0", h);
volVectorField hU0("hU0", hU);

volScalarField k1_h
(
    IOobject("k1_h", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField k2_h
(
    IOobject("k2_h", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField k3_h
(
    IOobject("k3_h", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField k4_h
(
    IOobject("k4_h", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

volVectorField k1_hU
(
    IOobject("k1_hU", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
);
volVectorField k2_hU
(
    IOobject("k2_hU", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
);
volVectorField k3_hU
(
    IOobject("k3_hU", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
);
volVectorField k4_hU
(
    IOobject("k4_hU", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
);
