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
    checkMesh

Description
    Checks validity of a mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
#include "foamTime.H"

#include "polyMesh.H"
#include "globalMeshData.H"

#include "printMeshStats.H"
#include "checkTopology.H"
#include "checkGeometry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
    argList::validOptions.insert("noTopology", "");
    argList::validOptions.insert("allGeometry", "");
    argList::validOptions.insert("allTopology", "");

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedPolyMesh.H"

    const bool noTopology  = args.optionFound("noTopology");
    const bool allGeometry = args.optionFound("allGeometry");
    const bool allTopology = args.optionFound("allTopology");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        polyMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            !timeI
         || state == polyMesh::TOPO_CHANGE
         || state == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Clear mesh before checking
            mesh.clearOut();

            // Reconstruct globalMeshData
            mesh.globalData();

            printMeshStats(mesh, allTopology);

            label noFailedChecks = 0;

            if (!noTopology)
            {
                noFailedChecks += checkTopology(mesh, allTopology, allGeometry);
            }

            noFailedChecks += checkGeometry(mesh, allGeometry);

            // Note: no reduction in noFailedChecks necessary since it is
            //       a counter of checks, not counter of failed cells,faces etc.

            if (noFailedChecks == 0)
            {
                Info<< nl << "Mesh OK." << nl << endl;
            }
            else
            {
                Info<< "\nFailed " << noFailedChecks << " mesh checks."
                    << nl << endl;
            }
        }
        else if (state == polyMesh::POINTS_MOVED)
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            label nFailedChecks = checkGeometry(mesh, allGeometry);

            if (nFailedChecks)
            {
                Info<< "\nFailed " << nFailedChecks << " mesh checks.\n"
                    << endl;
            }
            else
            {
                Info<< nl << "Mesh OK." << nl << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
