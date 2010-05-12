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
    checkMesh

Description
    Checks validity of a mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "globalMeshData.H"

#include "printMeshStats.H"
#include "checkTopology.H"
#include "checkGeometry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

#   include "addTimeOptionsNoConstant.H"

    argList::validOptions.insert("noTopology", "");
    argList::validOptions.insert("allGeometry", "");
    argList::validOptions.insert("allTopology", "");

#   include "setRootCase.H"

    const bool noTopology = args.options().found("noTopology");
    const bool allGeometry = args.options().found("allGeometry");
    const bool allTopology = args.options().found("allTopology");

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

#   include "checkTimeOptionsNoConstant.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createNamedPolyMesh.H"

    bool firstCheck = true;

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        polyMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            firstCheck
         || state == polyMesh::TOPO_CHANGE
         || state == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            firstCheck = false;

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

            reduce(noFailedChecks, sumOp<label>());

            if (noFailedChecks == 0)
            {
                Info<< "\nMesh OK."
                    << nl << endl;
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

            label noFailedChecks = checkGeometry(mesh, allGeometry);

            reduce(noFailedChecks, sumOp<label>());

            if (noFailedChecks == 0)
            {
                Info << "\nMesh OK."
                    << nl << endl;
            }
            else
            {
                Info<< "\nFailed " << noFailedChecks << " mesh checks."
                    << nl << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
