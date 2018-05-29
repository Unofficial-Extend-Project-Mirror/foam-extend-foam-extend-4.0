/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    calcOverset

Description
    Calculates overset addressing and everything necessary for overset
    interpolation. Used for parallel scaling tests of overset assembly.

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "oversetMesh.H"
#include "oversetFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    volVectorField cellCentres
    (
        IOobject
        (
            "cellCentres",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C(),
        "zeroGradient"
    );

    volVectorField origCellCentres("origCellCentres", cellCentres);
    origCellCentres.write();

    // Create CPU time object
    cpuTime oversetTime;

    // Make sure everything necessary for overset interpolation has been
    // calculated by calling a single overset interpolation
    oversetFvPatchVectorField::oversetInterpolate(cellCentres);

    Info<< "Elapsed time for overset assembly and single interpolation: "
        << oversetTime.cpuTimeIncrement() << endl;

    // Write interpolated cell centres
    cellCentres.write();

#   include "writeOversetMasks.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
