/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Update the polyMesh corresponding to the given map.

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "foamTime.H"
#include "globalMeshData.H"
#include "meshObjectBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Update zones.  Since boundary depends on zones, they need to be
    // updated first.  HJ, 20/May/2014
    pointZones_.updateMesh();
    faceZones_.updateMesh();
    cellZones_.updateMesh();

    // Update boundaryMesh (note that patches themselves are already ok)
    boundary_.updateMesh();

    // Clear out parallel data.  HJ, 27/Nov/2009
    deleteDemandDrivenData(globalMeshDataPtr_);

    setInstance(time().timeName());

    // Map the old motion points if present
    if (oldAllPointsPtr_)
    {
        // Make a copy of the original points
        pointField oldMotionPoints = *oldAllPointsPtr_;

        pointField& newMotionPoints = *oldAllPointsPtr_;

        // Resize the list to new size
        newMotionPoints.setSize(allPoints_.size());

        // Map the list
        newMotionPoints.map(oldMotionPoints, mpm.pointMap());

        // Reset old points if present
        if (oldPointsPtr_)
        {
            oldPointsPtr_->reset(*oldAllPointsPtr_, nPoints());
        }
    }

    // Reset valid directions (could change by faces put into empty patches)
    geometricD_ = Vector<label>::zero;
    solutionD_ = Vector<label>::zero;

    // Update all function objects
    // Moved from fvMesh.C in 1.6.x merge.  HJ, 29/Aug/2010
    meshObjectBase::allUpdateTopology<polyMesh>(*this, mpm);
}


// Sync mesh update with changes on other processors
void Foam::polyMesh::syncUpdateMesh()
{
    // Update zones.  Since boundary depends on zones, they need to be
    // updated first.  HJ, 20/May/2014
    pointZones_.updateMesh();
    faceZones_.updateMesh();
    cellZones_.updateMesh();

    // Update boundaryMesh (note that patches themselves already ok)
    boundary_.updateMesh();

    // Clear out parallel data.  HJ, 27/Nov/2009
    deleteDemandDrivenData(globalMeshDataPtr_);

    setInstance(time().timeName());

    // Reset valid directions (could change by faces put into empty patches)
    geometricD_ = Vector<label>::zero;
    solutionD_ = Vector<label>::zero;

    // Update all function objects
    // Moved from fvMesh.C in 1.6.x merge.  HJ, 29/Aug/2010

    // Instantiate a dummy mapPolyMesh
    autoPtr<mapPolyMesh> mapPtr(new mapPolyMesh(*this));

    meshObjectBase::allUpdateTopology<polyMesh>(*this, mapPtr());
}


// ************************************************************************* //
