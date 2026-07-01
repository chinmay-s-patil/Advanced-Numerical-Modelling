/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    shallowWaterFoam

Description
    Explicit solver for the linearized 2D Shallow Water Equations.

    Prognostic variables:
      zeta  -- free-surface height perturbation [m]
      qFlux -- depth-integrated transport (U, V) [m^2/s]

    Time integration: Euler, RK4, or Crank-Nicolson (predictor-corrector).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    #include "createRK4Fields.H"
    // #include "createCNFields.H"

    while (runTime.run())
    {
        #include "setDeltaT.H"

        runTime ++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (type == "euler")
        {
            #include "eulerStep.H"
        }

        else if (type == "rk4")
        {
            #include "rk4Step.H"
        }

        // else if (type == "cn")
        // {
        //     Info << "Explicit Crank Nicolson method selected" << endl;
        //     #include "cnStep.H"
        // }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
