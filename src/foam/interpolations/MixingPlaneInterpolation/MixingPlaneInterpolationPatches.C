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
    Mixing plane class dealing with transfer of data between two
    primitivePatches

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "MixingPlaneInterpolationTemplate.H"
#include "demandDrivenData.H"
#include "PrimitivePatchTemplate.H"
#include "IOmanip.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Find and fix the patch faces that have cylindrical points coordinates in
//  contact with or split on both sides of the -180/+180 degree
//  axis eg: [-179,+179]
//
//  The cure is simple: we first find the problematic faces, (there might be
//  none), then we modify only their angle coordinate in
//  order to bring the whole face across the angle = +0 axis instead.
//  This is a simple shift in cylindrical coordinates,
//  but also a pure rotation in cartesian space for the faces.
//
//  Since the only objective for these localCoordFaces is to create a patch for
//  constructing a temporary GGI in cylindrical coordinate
//  space, a single shift will have no effect whatsoever on the GGI weights
//   because we are using 360 degrees face ribbons for the other
//  GGI patch.
//
//  By doing so, we will obvioulsy mess up the clean topology of the original
//  cylindrical coordinates patch by shifting around a few faces;
//  but keep in mind the only topology necessary for the GGI is the faces
//  topology. A bunch of arbitrary located faces is as good as a bunch
//  of cleanly laid-out faces into a regular patch. By only shifting the face
//  along the angle axis, the face topology is preserved.
//
// Another interesting side-effect of this procedure is that it will cut
//  a 360 degrees patch along the -180/+180 axis, leaving two un-connected
// "right-side" and "left-side" edges when switching the patch points to
//  cylindrical coordinates.
template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::correctStraddlingFaces
(
    faceList& localCoordFaces,
    pointField& localCoordFacesPoint
) const
{
    // Memorize the list of displaced points so we can add them once at the
    // end. This minimizes the
    // memory reallocation necessary if we would do it one at a time instead.

    typedef std::map<label, point> labelPointMap;

    labelPointMap newFacePoints;

    label newPointLabel = localCoordFacesPoint.size();

    // Collect data on span limiting from coordinate system
    coordinateSystem::spanInfo spanLimited = cs_.spanLimited();
    boundBox spanBounds = cs_.spanBounds();

    const direction sweepDir = sweepAxisSwitch();

    // If there are no limits on both sides of the span,
    // adjustment can be skipped
    if
    (
        !spanLimited[sweepDir].first()
     || !spanLimited[sweepDir].second()
    )
    {
        // Nothing to do
        return;
    }

    // Check if coordinate system is limited in spanwise range
    scalar maxSpan = spanBounds.max().component(sweepDir);

    if (debug)
    {
        Info<< "Fixing span straddle in direction " << sweepDir
            << " for span " << maxSpan << endl;
    }

    // Do all faces
    forAll (localCoordFaces, sFi)
    {
        bool faceIsOk(true);

        scalarField pointsAngleValues =
            localCoordFaces[sFi].points
            (
                localCoordFacesPoint
            ).component(sweepDir);

        for (label i = 0; i < pointsAngleValues.size(); i++)
        {
            for (label j = i + 1; j < pointsAngleValues.size(); j++)
            {
                // We shift away any faces in contact with the
                // straddling axis
                if
                (
                    (
                        mag(pointsAngleValues[i] - pointsAngleValues[j])
                      > maxSpan
                    )
                 || (mag(mag(pointsAngleValues[i]) - maxSpan)) < VSMALL
                 || (mag(mag(pointsAngleValues[j]) - maxSpan)) < VSMALL
                )
                {
                    // We need to correct this
                    // Grab the original points labels
                    labelList pointLbls = localCoordFaces[sFi];

                    forAll (pointLbls, ptI)
                    {
                        point ptCoord = localCoordFacesPoint[pointLbls[ptI]];

                        // Switch point across the +0 angle axis
                        if (ptCoord[sweepDir] < 0.0)
                        {
                            ptCoord[sweepDir] = ptCoord[sweepDir] + maxSpan;
                        }
                        else
                        {
                            ptCoord[sweepDir] = ptCoord[sweepDir] - maxSpan;
                        }

                        // Memorize data in order to reallocate
                        // localCoordFacesPoint only once
                        newFacePoints.insert
                        (
                            std::pair<label, point>(newPointLabel, ptCoord)
                        );

                        // Memorize the new point label right away for this
                        // face and increment index for the next
                        localCoordFaces[sFi][ptI] = newPointLabel++;
                    }
                    // Done with this face
                    faceIsOk = false;

                    break;
                }
            }

            if (!faceIsOk) break;
        }
    }

    // Transfer the new points information into localCoordFacesPoint using
    // only one reallocation
    // We probably have point duplicates here; no big deal...
    localCoordFacesPoint.setSize
    (
        localCoordFacesPoint.size() + newFacePoints.size()
    );

    forAllIter(labelPointMap, newFacePoints, nPi)
    {
          // First == label, Second == point
        localCoordFacesPoint[nPi->first] = nPi->second;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
calcTransformedPatches() const
{
    if
    (
        transformedMasterPatchPtr_
     || transformedShadowPatchPtr_
    )
    {
        FatalErrorIn
        (
            "void MixingPlaneInterpolation::calcTransformedPatches() const"
        )   << "Patches already calculated"
            << abort(FatalError);
    }

    // Duplicate the master/slave patch faces
    faceList masterFaces = masterPatch_.localFaces();
    faceList slaveFaces = slavePatch_.localFaces();

    // Let's compute the patches coordinates values in mixing coordinates
    pointField masterPointsLocalCoord =
        cs_.localPosition(masterPatch_.localPoints());

    pointField slavePointsLocalCoord =
        cs_.localPosition(slavePatch_.localPoints());

    // Next, we need to find and fix the patch faces that have straddled
    // the span
    // Note: The face indices and the local coordinates and the
    // number of points may be modified within correctStraddlingFaces
    // MB, 27/Jan/2011

    correctStraddlingFaces(masterFaces, masterPointsLocalCoord);
    correctStraddlingFaces(slaveFaces, slavePointsLocalCoord);

    if(debug)
    {
        InfoIn
        (
            "MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "calcTransformedPatches()"
        )   << "masterPointsLocalCoord: "
            << masterPointsLocalCoord << nl
            << "slavePointsLocalCoord: "
            << slavePointsLocalCoord << endl;
    }

    // Create the local coords patches

    transformedMasterPatchPtr_ =
        new standAlonePatch
        (
            masterFaces,
            masterPointsLocalCoord
        );

    transformedShadowPatchPtr_ =
        new standAlonePatch
        (
            slaveFaces,
            slavePointsLocalCoord
        );
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::calcMixingPlanePatch() const
{
    if (mixingPlanePatchPtr_)
    {
        FatalErrorIn
        (
            "void MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "calcMixingPlanePatch() const"
        )   << "Circumferential average patch alreacy calculated"
            << abort(FatalError);
    }

    const pointField& profile = this->interpolationProfile();

    // Create patch from profile

    // Grab the bounding boxes of master and slave

    boundBox masterBB
    (
        transformedMasterPatch().localPoints(),
        false
    );

    boundBox slaveBB
    (
        transformedShadowPatch().localPoints(),
        false
    );

    // Collect data on span limiting from coordinate system
    coordinateSystem::spanInfo spanLimited = cs_.spanLimited();
    boundBox spanBounds = cs_.spanBounds();

    const direction sweepDir = sweepAxisSwitch();

    // Get span from bounding boxes
    scalar minSpan =
        Foam::min
        (
            masterBB.min()[sweepDir],
            slaveBB.min()[sweepDir]
        ) - SMALL;

    scalar maxSpan =
        Foam::max
        (
            masterBB.max()[sweepDir],
            slaveBB.max()[sweepDir]
        ) + SMALL;

    if (debug)
    {
        InfoIn
        (
            "MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "calcMixingPlanePatch() const"
        )   << "minSpan from patch BB : " << minSpan << nl
            << "maxSpan from patch BB : " << maxSpan << endl;
    }

    // Correct for limited span
    if (spanLimited[sweepDir].first())
    {
        minSpan = spanBounds.min()[sweepDir];
    }

    if (spanLimited[sweepDir].second())
    {
        maxSpan = spanBounds.max()[sweepDir];
    }

    if (debug)
    {
        InfoIn
        (
            "MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "calcMixingPlanePatch() const"
        )   << "minSpan after checking spanLimited: " << minSpan << nl
            << "maxSpan after checking spanLimited: " << maxSpan << endl;
    }

    label nRibbons = profile.size() - 1;

    // Points for averaging patch in the local coordinate system
    pointField mixingPatchPoints(2*nRibbons + 2);
    label nextPointID = 0;

    faceList mixingPatchFaces(nRibbons);
    scalarField mixingPatchFacesArea(nRibbons);

    // Insert lower bound points and expand to bounds
    mixingPatchPoints[0] = profile[0];
    mixingPatchPoints[0].replace(sweepDir, minSpan);

    mixingPatchPoints[1] = profile[0];
    mixingPatchPoints[1].replace(sweepDir, maxSpan);
    nextPointID = 2;


    // Create patch faces
    forAll (mixingPatchFaces, fI)
    {
        label fi2 = 2*fI;
        // Add top bound points and expand to bounds
        mixingPatchPoints[nextPointID] = profile[fI + 1];
        mixingPatchPoints[nextPointID].replace(sweepDir, minSpan);
        nextPointID++;

        mixingPatchPoints[nextPointID] = profile[fI + 1];
        mixingPatchPoints[nextPointID].replace(sweepDir, maxSpan);
        nextPointID++;

        face curFace(4);
        curFace[0] = fi2;
        curFace[1] = fi2 + 1;
        curFace[2] = fi2 + 3;
        curFace[3] = fi2 + 2;

        mixingPatchFaces[fI] = curFace;
    }

    mixingPlanePatchPtr_ =
        new standAlonePatch
        (
            mixingPatchFaces,
            mixingPatchPoints
        );


    if (debug > 0)
    {
        InfoIn
        (
            "MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "calcMixingPlanePatch() const"
        )   << "mixingPatch: "
            << *mixingPlanePatchPtr_ << nl
            << "mixingPatch.points : "
            << mixingPlanePatchPtr_->points() << endl;
    }
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::clearTransfomedPatches()
{
    deleteDemandDrivenData(transformedMasterPatchPtr_);
    deleteDemandDrivenData(transformedShadowPatchPtr_);
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::clearMixingPlanePatch()
{
    deleteDemandDrivenData(mixingPlanePatchPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
