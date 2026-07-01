// setDeltaT.H -- CFL-based time step for linearized SWE
// Wave speed: c = sqrt(g * H)
// Only adjusts deltaT if adjustTimeStep is enabled in controlDict.
{
    // Wave speed at faces
    surfaceScalarField cSf
    (
        fvc::interpolate
        (
            sqrt(magg * max(H, dimensionedScalar("Hmin", dimLength, 1e-6)))
        ) * mesh.magSf()
    );

    scalarField sumCSf
    (
        fvc::surfaceSum(cSf)().primitiveField()
    );

    scalar CoNum =
        0.5*gMax(sumCSf/mesh.V().field())*runTime.deltaTValue();

    Info<< "Max wave speed Courant Number = " << CoNum
        << "  deltaT = " << runTime.deltaTValue() << endl;

    if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        scalar deltaT = min
        (
            maxCo*runTime.deltaTValue()/max(CoNum, SMALL),
            maxDeltaT
        );

        runTime.setDeltaT(deltaT);
    }
}
