{
    surfaceScalarField amaxSf
    (
        fvc::interpolate(mag(U) + sqrt(magg*max(h, dimensionedScalar("hMin", dimLength, 1e-6))))*mesh.magSf()
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
