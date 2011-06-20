/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    multiSolverDemo

Description
    Combination of icoFoam and scalarTransportFoam for testing of new multiTime
    framework.

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createMultiSolver.H"

// * * * * * * * * * * * * * * * * icoFoam1  * * * * * * * * * * * * * * * * //

    Info << "*** Switching to icoFoam1 ***\n" << endl;
    solverDomain = "icoFoam1";
#   include "setSolverDomain.H"

#   include "solverIcoFoam.H"

// * * * * * * * * * * * * * scalarTransportFoam * * * * * * * * * * * * * * //

    Info << "*** Switching to scalarTransportFoam ***\n" << endl;
    solverDomain = "scalarTransportFoam";
#   include "setSolverDomain.H"

#   include "solverScalarTransportFoam.H"

    multiRun++;
    
// * * * * * * * * * * * * * * * * icoFoam2  * * * * * * * * * * * * * * * * //

    Info << "*** Switching to icoFoam2 ***\n" << endl;
    solverDomain = "icoFoam2";
#   include "setSolverDomain.H"

#   include "solverIcoFoam.H"

// * * * * * * * * * * * * * scalarTransportFoam * * * * * * * * * * * * * * //

    Info << "*** Switching to scalarTransportFoam ***\n" << endl;
    solverDomain = "scalarTransportFoam";
#   include "setSolverDomain.H"

#   include "solverScalarTransportFoam.H"

#   include "endMultiSolver.H"
    return(0);
}

// ************************************************************************* //
