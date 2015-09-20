/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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

    // Only necessary if we revisit the same solver domain twice in the same
    // superLoop (scalarTransportFoam, in this case)
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
