// Crank Nicolson Method

volScalarField hStar("hStar", h);
volScalarField RHS_h("RHS_h", h*0);
volScalarField RHS_h_Star("RHS_h_Star", h*0);

volVectorField UStar("UStar", U);

volVectorField hUStar("hUStar", hU);
volVectorField RHS_hU("RHS_hU", hU*0);
volVectorField RHS_hU_Star("RHS_hU_Star", hU*0);
