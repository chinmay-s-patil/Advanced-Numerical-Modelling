const dictionary& timeControls = runTime.controlDict();

scalar maxCo
(
    timeControls.lookupOrDefault<scalar>("maxCo", 0.5)
);

scalar maxDeltaT
(
    timeControls.lookupOrDefault<scalar>("maxDeltaT", GREAT)
);
