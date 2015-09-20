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
    boxTurb3D

Description
    Makes a box of turbulence which conforms to a given energy
    spectrum and is divergence free.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "graph.H"
#include "OFstream.H"
#include "Kmesh.H"
#include "turbGen.H"
#include "calcEk.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readBoxTurbDict.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Kmesh K(mesh);

    turbGen Ugen(K, Ea, k0);

    U.internalField() = Ugen.U();
    U.correctBoundaryConditions();

    Info<< "k("
         << runTime.timeName()
         << ") = "
         << 3.0/2.0*average(magSqr(U)).value() << endl;

    U.write();

    calcEk(U, K).write(runTime.timePath()/"Ek", runTime.graphFormat());

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
