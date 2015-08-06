/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenChecks.H"
#include "polyMeshGenAddressing.H"
#include "pyramidPointFaceRef.H"
#include "tetrahedron.H"
#include "cell.H"
#include "mathematicalConstants.H"
#include "ListOps.H"
#include "Map.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace polyMeshGenChecks
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool checkGeometry(const polyMeshGen& mesh, const bool report)
{
    label noFailedChecks(0);

    if( checkClosedBoundary(mesh, report) ) ++noFailedChecks;
    if( checkClosedCells(mesh, report) ) ++noFailedChecks;
    if( checkFaceAreas(mesh, report) ) ++noFailedChecks;
    if( checkCellVolumes(mesh, report) ) ++noFailedChecks;
    if( checkFaceDotProduct(mesh, report) ) ++noFailedChecks;
    if( checkFaceUniformity(mesh, report) ) ++noFailedChecks;
    if( checkFacePyramids(mesh, report) ) ++noFailedChecks;
    if( checkFaceSkewness(mesh, report) ) ++noFailedChecks;
    if( checkCellPartTetrahedra(mesh, report) ) ++noFailedChecks;

    if( noFailedChecks == 0 )
    {
        if( report )
            Info << "Mesh geometry OK." << endl;

        return false;
    }
    else
    {
        Info<< "Failed " << noFailedChecks << " mesh geometry checks." << endl;

        return true;
    }
}

bool checkTopology(const polyMeshGen& mesh, const bool report)
{
    label noFailedChecks(0);

    if( checkPoints(mesh, report) ) ++noFailedChecks;
    if( checkUpperTriangular(mesh, report) ) ++noFailedChecks;
    if( checkCellsZipUp(mesh, report) ) ++noFailedChecks;
    if( checkFaceVertices(mesh, report) ) ++noFailedChecks;

    if( noFailedChecks == 0 )
    {
        if( report )
            Info << "Mesh topology OK." << endl;

        return false;
    }
    else
    {
        Info<< "Failed " << noFailedChecks << " mesh topology checks." << endl;

        return true;
    }
}

bool checkMesh(const polyMeshGen& mesh, const bool report)
{
    bool failedChecks = checkTopology(mesh, report);
    failedChecks |= checkGeometry(mesh, report);

    if( !failedChecks )
    {
        if( report )
            Info << "Mesh OK." << endl;

        return false;
    }
    else
    {
        Info << "Failed some mesh checks." << endl;

        return true;
    }
}

} // End namespace polyMeshGenChecks

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
