{
    surfaceScalarField amaxSf
    (
        fvc::interpolate(mag(U) + sqrt(magg*h))*mesh.magSf()
    );

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    scalar CoNum =
        0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();

    scalar deltaT = min
    (
        maxCo*runTime.deltaTValue()/max(CoNum, SMALL),
        maxDeltaT
    );

    runTime.setDeltaT(deltaT);

    Info<< "Max wave speed Courant Number = " << CoNum
        << "  deltaT = " << runTime.deltaTValue() << endl;
}
