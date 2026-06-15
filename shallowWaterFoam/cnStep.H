// Crank Nicolson Method
// 
// 1. Compute RHS at current state (n) -> RHS_h, RHS_hU
// 2. Euler step -> hStar = h + dt*RHS_h, hUStar = hU + dt*RHS_hU
// 3. Compute RHS at star state -> RHS_h_star, RHS_hU_star
// 4. Average -> h = h + dt/2*(RHS_h + RHS_h_star)
//              hU = hU + dt/2*(RHS_hU + RHS_hU_star)
// 5. Recover U = hU/h

surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(hU) & mesh.Sf()
);

// Snapshoting h and hU
RHS_h = fvc::div(-phi);
RHS_hU = -fvc::div(phi, U)
       - (magg/2)*fvc::grad(h*h)
       - magg*h*fvc::grad(zb);

// Step 2: Euler Predictor
hStar = h + runTime.deltaT()*RHS_h;
hUStar = hU + runTime.deltaT()*RHS_hU;
UStar = hUStar/hStar;

// Step 3: RHS at star state
surfaceScalarField phiStar
(
    IOobject("phiStar", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(hUStar) & mesh.Sf()
);

// Computing Star Values for h and hU
RHS_h_Star = fvc::div(-phiStar);
RHS_hU_Star =   -fvc::div(phiStar, UStar)
                - (magg/2)*fvc::grad(hStar*hStar)
                - magg*hStar*fvc::grad(zb);


// Step 4: Averaging Values
h = h + (runTime.deltaT()/2)*(RHS_h + RHS_h_Star);
hU = hU + (runTime.deltaT()/2)*(RHS_hU + RHS_hU_Star);

// Step 5: Recover U
U = hU/h;
