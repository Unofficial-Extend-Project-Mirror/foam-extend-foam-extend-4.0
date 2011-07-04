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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "overlapGgiPolyPatch.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "OFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Create expanded patch
Foam::standAlonePatch*
Foam::overlapGgiPolyPatch::calcExpandedGeometry(label ncp, label index) const
{
    const scalar myAngle = 360.0/scalar(ncp);

    // Create expanded master points and faces
    const faceZone& geomZone = boundaryMesh().mesh().faceZones()[index];
    const primitiveFacePatch& geom = geomZone();

    const pointField& geomLocalPoints = geom.localPoints();

    pointField expandedPoints(ncp*geomLocalPoints.size());

    // Transform points
    label nPointsGeom = 0;

    for (label copyI = 0; copyI < ncp; copyI++)
    {
        // Calculate transform
        const tensor curRotation =
            RodriguesRotation(rotationAxis_,  copyI*myAngle);

        forAll (geomLocalPoints, pointI)
        {
            expandedPoints[nPointsGeom] =
                transform(curRotation, geomLocalPoints[pointI]);
            nPointsGeom++;
        }
    }

    // Transform faces
    const faceList& geomLocalFaces = geom.localFaces();
    faceList expandedFaces(ncp*geomLocalFaces.size());

    label nFacesGeom = 0;

    for (label copyI = 0; copyI < ncp; copyI++)
    {
        const label copyOffsetGeom = copyI*geomLocalPoints.size();

        forAll (geomLocalFaces, faceI)
        {
            const face& curGeomFace = geomLocalFaces[faceI];

            face& curExpandedFace = expandedFaces[nFacesGeom];

            // Copy face with offsets
            curExpandedFace.setSize(curGeomFace.size());

            forAll (curGeomFace, fpI)
            {
                curExpandedFace[fpI] = curGeomFace[fpI] + copyOffsetGeom;
            }

            nFacesGeom++;
        }
    }

    if (debug > 1)
    {
        Info << "Writing expanded geom patch as VTK" << endl;

        const polyMesh& mesh = boundaryMesh().mesh();

        fileName fvPath(mesh.time().path()/"VTK");
        mkDir(fvPath);

        OStringStream outputFilename;
        outputFilename << "expandedGeom" << name() << shadow().name()
            << index << "_" << mesh.time().timeName();

        standAlonePatch::writeVTK
        (
            fvPath/fileName(outputFilename.str()),
            expandedFaces,
            expandedPoints
        );
    }

    return new standAlonePatch(expandedFaces, expandedPoints);
}


const Foam::standAlonePatch& Foam::overlapGgiPolyPatch::expandedMaster() const
{
    if (!expandedMasterPtr_)
    {
        expandedMasterPtr_ = calcExpandedGeometry(nCopies(), zoneIndex());
    }

    return *expandedMasterPtr_;
}


const Foam::standAlonePatch& Foam::overlapGgiPolyPatch::expandedSlave() const
{
    if (!expandedSlavePtr_)
    {
        expandedSlavePtr_ =
            calcExpandedGeometry(shadow().nCopies(), shadow().zoneIndex());
    }

    return *expandedSlavePtr_;
}


void Foam::overlapGgiPolyPatch::calcPatchToPatch() const
{
    // Create patch-to-patch interpolation between the expanded master
    // and slave patches
    if (patchToPatchPtr_)
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcPatchToPatch() const")
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        patchToPatchPtr_ =
            new overlapGgiInterpolation
            (
                expandedMaster(),
                expandedSlave(),
                forwardT(),
                reverseT(),
                separation(),
                0,             // master overlap tolerance
                0,             // slave overlap tolerance
                true,          // Rescale weighting factors.  Bug fix, MB.
//                 ggiInterpolation::AABB
                overlapGgiInterpolation::BB_OCTREE  // Octree search, MB.

            );

        // Abort immediatly if uncovered faces are present
        if
        (
            patchToPatch().uncoveredMasterFaces().size() > 0
         || patchToPatch().uncoveredSlaveFaces().size() > 0
        )
        {
            FatalErrorIn("void overlapGgiPolyPatch::calcPatchToPatch() const")
                << "Found uncovered faces for GGI interface "
                << name() << "/" << shadowName() << endl
                << "This is an unrecoverable error. Aborting."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcPatchToPatch() const")
            << "Attempting to create GGIInterpolation on a shadow"
            << abort(FatalError);
    }
}


const Foam::overlapGgiInterpolation&
Foam::overlapGgiPolyPatch::patchToPatch() const
{
    if (master())
    {
        if (!patchToPatchPtr_)
        {
            if (debug)
            {
                Info<< "Initializing the GGI interpolator between "
                    << "master/shadow patches: "
                    << name() << "/" << shadowName()
                    << endl;
            }

            calcPatchToPatch();
        }

        return *patchToPatchPtr_;
    }
    else
    {
        return shadow().patchToPatch();
    }
}


void Foam::overlapGgiPolyPatch::calcReconFaceCellCentres() const
{
    if (reconFaceCellCentresPtr_)
    {
        FatalErrorIn
        (
            "void overlapGgiPolyPatch::calcReconFaceCellCentres() const"
        )   << "Reconstructed cell centres already calculated"
            << abort(FatalError);
    }

    // Create neighbouring face centres using interpolation
    if (master())
    {
        const label shadowID = shadowIndex();

        // Get the transformed and interpolated shadow face cell centers
        reconFaceCellCentresPtr_ =
            new vectorField
            (
                interpolate
                (
                    boundaryMesh()[shadowID].faceCellCentres()
                  - boundaryMesh()[shadowID].faceCentres()
                )
              + faceCentres()
            );
    }
    else
    {
        FatalErrorIn
        (
            "void overlapGgiPolyPatch::calcReconFaceCellCentres() const"
        )   << "Attempting to create reconFaceCellCentres on a shadow"
            << abort(FatalError);
    }
}


void Foam::overlapGgiPolyPatch::checkDefinition() const
{
    // Sanity checks
    // 1. Check
    Info << "overlapGgiPolyPatch: sanity checks missing.  HJ" << endl;

//     if
//     (
//         (rotationAngle() - shadow().rotationAngle()) > SMALL
//      || cmptSum(rotationAxis() - shadow().rotationAxis()) > SMALL
//     )
//     {
//         FatalErrorIn("overlapGgiPolyPatch::patchToPatch")
//             << "    The rotation angle for patch name           : "
//             << name() << " is: " << rotationAngle() << " axis: "
//             << shadow().rotationAxis() << nl
//             << "    The rotation angle for the shadow patch name: "
//             << shadowName() << " is: "
//             << shadow().rotationAngle() << " axis: "
//             << shadow().rotationAxis() << nl
//             << "    Both values need to be the same in "
//             << "the boundary file "
//             << abort(FatalError);
//     }
}


void Foam::overlapGgiPolyPatch::clearGeom()
{
    // Note: Currently deleting patch-to-patch interpolation together with
    // expanded master and slave patches on mesh motion to avoid problems
    // with motion of points in primitive patch.
    // HJ, 4/Jul/2011
    deleteDemandDrivenData(expandedMasterPtr_);
    deleteDemandDrivenData(expandedSlavePtr_);

    deleteDemandDrivenData(patchToPatchPtr_);

    deleteDemandDrivenData(reconFaceCellCentresPtr_);
}


void Foam::overlapGgiPolyPatch::clearOut()
{
    clearGeom();

    deleteDemandDrivenData(localParallelPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField&
Foam::overlapGgiPolyPatch::reconFaceCellCentres() const
{
    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }

    return *reconFaceCellCentresPtr_;
}


void Foam::overlapGgiPolyPatch::initAddressing()
{
    if (active())
    {
        // Calculate transforms for correct GGI interpolator cut
        calcTransforms();
        localParallel();
    }

    polyPatch::initAddressing();
}


void Foam::overlapGgiPolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::overlapGgiPolyPatch::initGeometry()
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011
    if (active())
    {
        // Note: Only master calculates recon; slave uses master interpolation
        if (master())
        {
            reconFaceCellCentres();
        }
    }

    polyPatch::initGeometry();
}


void Foam::overlapGgiPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}


void Foam::overlapGgiPolyPatch::initMovePoints(const pointField& p)
{
    clearGeom();

    // Calculate transforms on mesh motion?
    calcTransforms();

    // Update interpolation for new relative position of GGI interfaces
    // Note: currently, patches and interpolation are cleared in clearGeom()
    // HJ. 4/Jul/2011
//     if (patchToPatchPtr_)
//     {
//         patchToPatchPtr_->movePoints();
//     }

    if (active() && master())
    {
        reconFaceCellCentres();
    }

    polyPatch::initMovePoints(p);
}


void Foam::overlapGgiPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
}


void Foam::overlapGgiPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::overlapGgiPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    clearOut();
}


void Foam::overlapGgiPolyPatch::calcTransforms()
{
    forwardT_.setSize(0);
    reverseT_.setSize(0);
    separation_.setSize(0);
}


void Foam::overlapGgiPolyPatch::initOrder(const primitivePatch& pp) const
{}


bool Foam::overlapGgiPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size(), -1);
    rotation.setSize(pp.size(), 0);

    // Nothing changes
    return false;
}


void Foam::overlapGgiPolyPatch::syncOrder() const
{}


// ************************************************************************* //
