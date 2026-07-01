// createRK4Fields.H -- RK4 intermediate storage for linearized SWE
//
// zeta has dimensions [m],        d/dt(zeta) has dimensions [m/s]
// qFlux has dimensions [m^2/s],   d/dt(qFlux) has dimensions [m^2/s^2]

volScalarField  zeta0("zeta0", zeta);
volVectorField  qFlux0("qFlux0", qFlux);

// k-stages for zeta (dzeta/dt ~ [m/s])
volScalarField k1_zeta
(
    IOobject("k1_zeta", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField k2_zeta("k2_zeta", k1_zeta);
volScalarField k3_zeta("k3_zeta", k1_zeta);
volScalarField k4_zeta("k4_zeta", k1_zeta);

// k-stages for qFlux (dQ/dt ~ [m^2/s^2])
volVectorField k1_qFlux
(
    IOobject("k1_qFlux", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
);
volVectorField k2_qFlux("k2_qFlux", k1_qFlux);
volVectorField k3_qFlux("k3_qFlux", k1_qFlux);
volVectorField k4_qFlux("k4_qFlux", k1_qFlux);
