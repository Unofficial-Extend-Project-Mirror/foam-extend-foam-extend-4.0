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

Description

\*---------------------------------------------------------------------------*/

#include "processorTetPolyPatchFaceDecomp.H"
#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "globalTetPolyPatchFaceDecomp.H"
#include "tetPolyMeshFaceDecomp.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorTetPolyPatchFaceDecomp::calcMeshPoints() const
{
    if (meshPointsPtr_)
    {
        FatalErrorIn
        (
            "void processorTetPolyPatchFaceDecomp::calcMeshPoints() const"
        )   << "meshPointsPtr_ already allocated"
            << abort(FatalError);
    }

    // Algorithm:
    // Depending on whether the patch is a master or a slave, get the primitive
    // patch points and filter away the points from the global patch.

    labelList mp(0);

    if (isMaster())
    {
        mp = procPolyPatch_.meshPoints();
    }
    else
    {
        // Slave side.  Create the reversed patch and pick up its points
        // so that the order is correct
        const polyPatch& pp = patch();

        faceList masterFaces(pp.size());

        forAll (pp, faceI)
        {
            masterFaces[faceI] = pp[faceI].reverseFace();
        }

        mp = primitiveFacePatch
            (
                masterFaces,
                pp.points()
            ).meshPoints();
    }

    // Get reference to shared processor points
    const labelList& sharedPoints =
        boundaryMesh().globalPatch().meshPoints();

    // Filter the shared points out of the list
    meshPointsPtr_ = new labelList(mp.size() + procPolyPatch_.size());
    labelList& filtPoints = *meshPointsPtr_;

    label noFiltPoints = 0;

    forAll (mp, pointI)
    {
        label curP = mp[pointI];

        bool found = false;
        forAll (sharedPoints, sharedI)
        {
            if (sharedPoints[sharedI] == curP)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            filtPoints[noFiltPoints] = curP;
            noFiltPoints++;
        }
    }

    // insert faces using the offset
    // These cannot be shared
    const label faceOffset = boundaryMesh().mesh().faceOffset();

    const label polyPatchStart = patch().start() + faceOffset;
    const label polyPatchStartPlusSize = polyPatchStart + patch().size();

    for (label i = polyPatchStart; i < polyPatchStartPlusSize; i++)
    {
        filtPoints[noFiltPoints] = i;
        noFiltPoints++;
    }

    filtPoints.setSize(noFiltPoints);
}


// ************************************************************************* //
