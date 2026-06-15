// RK4 Step

// Snapshoting h and hU
h0 = h;
hU0 = hU;

// Stage 1: k1 at current state
surfaceScalarField phi1(fvc::interpolate(hU) & mesh.Sf());

k1_h  = -fvc::div(phi1);
k1_hU = -fvc::div(phi1, U)
       - (magg/2)*fvc::grad(h*h)
       - magg*h*fvc::grad(zb);

// Stage 2: k2 at half-step
h  = h0 + (runTime.deltaT()/2)*k1_h;
hU = hU0 + (runTime.deltaT()/2)*k1_hU;
U  = hU/h;

surfaceScalarField phi2(fvc::interpolate(hU) & mesh.Sf());

k2_h  = -fvc::div(phi2);
k2_hU = -fvc::div(phi2, U)
       - (magg/2)*fvc::grad(h*h)
       - magg*h*fvc::grad(zb);

// Stage 3: k3 at half-step
h  = h0 + (runTime.deltaT()/2)*k2_h;
hU = hU0 + (runTime.deltaT()/2)*k2_hU;
U  = hU/h;

surfaceScalarField phi3(fvc::interpolate(hU) & mesh.Sf());

k3_h  = -fvc::div(phi3);
k3_hU = -fvc::div(phi3, U)
       - (magg/2)*fvc::grad(h*h)
       - magg*h*fvc::grad(zb);

// Stage 4: k4 at half-step
h  = h0 + (runTime.deltaT())*k3_h;
hU = hU0 + (runTime.deltaT())*k3_hU;
U  = hU/h;

surfaceScalarField phi4(fvc::interpolate(hU) & mesh.Sf());

k4_h  = -fvc::div(phi4);
k4_hU = -fvc::div(phi4, U)
       - (magg/2)*fvc::grad(h*h)
       - magg*h*fvc::grad(zb);

// Final Stage: Updating h and hU
h = h0 + (runTime.deltaT()/6)*(k1_h + 2*k2_h + 2*k3_h + k4_h);
hU = hU0 + (runTime.deltaT()/6)*(k1_hU + 2*k2_hU + 2*k3_hU + k4_hU);

// Recover U
U = hU/h;
