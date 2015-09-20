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
    moveDynamicMesh

Description
    Mesh motion and topological mesh changes utility for an engine geometry

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
#include "engineTopoChangerMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("checkFrequency", "int");

#   include "setRootCase.H"
#   include "createEngineTime.H"
#   include "createEngineDynamicMesh.H"

    // Read check frequency
    label checkFrequency = 1;
    args.optionReadIfPresent("checkFrequency", checkFrequency);

    fileName path = runTime.caseName();
    OFstream volFile(path + "/totVol.Cyl");

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        volFile << runTime.timeName() << "\t" << sum(mesh.V()).value() << endl;

        if (isDir(runTime.path()/"VTK"))
        {
            Info << "Clear VTK directory" << endl;
            rmDir(runTime.path()/"VTK");
        }

        mesh.update();

#       include "checkVolContinuity.H"
#       include "meshCourantNo.H"

        if
        (
            checkEngineMesh
         && (runTime.timeIndex() % checkFrequency == 0)
        )
        {
            mesh.checkMesh(true);
        }

        volFile << runTime.timeName() << tab << sum(mesh.V()).value() << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
