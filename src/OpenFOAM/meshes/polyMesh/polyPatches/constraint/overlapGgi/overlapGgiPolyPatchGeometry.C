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

void Foam::overlapGgiPolyPatch::calcExpandedMaster() const
{
    // Create expanded master patch interpolation
    if (expandedMasterPtr_)
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcExpandedMaster() const")
            << "Expanded master already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        // Create expanded master patch
        const label ncpm = nCopies();

        Info << "Number of master copies: " << ncpm << endl;

        // Create expanded master points and faces
        const polyPatch& master = boundaryMesh()[index()];
        const pointField& masterLocalPoints = master.localPoints();

        pointField MasterExpandedPoints(ncpm*masterLocalPoints.size());

        const scalar masterAngle = angle();

        Info << "Master Angle is: " << masterAngle << endl;

        // Transform points
        label nPoints_master = 0;

        for (label copyI = 0; copyI < ncpm; copyI++)
        {
            // Calculate transform
            const tensor curRotation =
                RodriguesRotation(rotationAxis_,  copyI*(masterAngle));

            forAll (masterLocalPoints, pointI)
            {
                MasterExpandedPoints[nPoints_master] =
                    transform(curRotation, masterLocalPoints[pointI]);
                nPoints_master++;
            }
        }

        // Transform faces
        const faceList& masterLocalFaces = master.localFaces();
        faceList MasterExpandedFaces(ncpm*masterLocalFaces.size());

        label nFacesMaster = 0;

        for (label copyI = 0; copyI < ncpm; copyI++)
        {
            const label copyOffsetMaster = copyI*masterLocalPoints.size();

            forAll (masterLocalFaces, faceI)
            {
                const face& curMasterFace = masterLocalFaces[faceI];

                face& MastercurExpandedFace =
                    MasterExpandedFaces[nFacesMaster];

                // Copy face with offsets
                MastercurExpandedFace.setSize(curMasterFace.size());

                forAll (curMasterFace, fpI)
                {
                    MastercurExpandedFace[fpI] =
                        curMasterFace[fpI] + copyOffsetMaster;
                }

                nFacesMaster++;
            }
        }

        expandedMasterPtr_ =
            new standAlonePatch(MasterExpandedFaces, MasterExpandedPoints);

        if (debug > 1)
        {
            Info << "Writing expanded master patch as VTK" << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            standAlonePatch::writeVTK
            (
                fvPath/fileName("expandedMaster" + name() + shadow().name()),
                MasterExpandedFaces,
                MasterExpandedPoints
            );
        }
    }
    else
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcExpandedMaster() const")
            << "Attempting to create expanded master on a shadow"
            << abort(FatalError);
    }
}

void Foam::overlapGgiPolyPatch::calcExpandedSlave() const
{
    // Create expanded slave patch interpolation
    if (expandedSlavePtr_)
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcExpandedSlave() const")
            << "Expanded slave already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        // Create expanded patch
        const label ncp = shadow().nCopies();

        Info << "Number of slave copies: " << ncp << endl;

        // Create expanded points and faces
        const polyPatch& slave = boundaryMesh()[shadowIndex()];
        const pointField& slaveLocalPoints = slave.localPoints();

        pointField expandedPoints(ncp*slaveLocalPoints.size());

        const scalar slaveAngle = shadow().angle();

        // Transform points
        label nPoints = 0;

        for (label copyI = 0; copyI < ncp; copyI++)
        {
            // Calculate transform
            const tensor curRotation =
                RodriguesRotation(rotationAxis_,  copyI*slaveAngle);

            forAll (slaveLocalPoints, pointI)
            {
                expandedPoints[nPoints] =
                    transform(curRotation, slaveLocalPoints[pointI]);
                nPoints++;
            }
        }

        // Transform faces

        const faceList& slaveLocalFaces = slave.localFaces();
        faceList expandedFaces(ncp*slaveLocalFaces.size());

        label nFaces = 0;

        for (label copyI = 0; copyI < ncp; copyI++)
        {
            const label copyOffset = copyI*slaveLocalPoints.size();

            forAll (slaveLocalFaces, faceI)
            {
                const face& curSlaveFace = slaveLocalFaces[faceI];

                face& curExpandedFace = expandedFaces[nFaces];

                // Copy face with offsets
                curExpandedFace.setSize(curSlaveFace.size());

                forAll (curSlaveFace, fpI)
                {
                    curExpandedFace[fpI] = curSlaveFace[fpI] + copyOffset;
                }

                nFaces++;
            }
        }

        expandedSlavePtr_ = new standAlonePatch(expandedFaces, expandedPoints);

        if (debug > 1)
        {
            Info << "Writing expanded slave patch as VTK" << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            standAlonePatch::writeVTK
            (
                fvPath/fileName("expandedSlave" + name() + shadow().name()),
                expandedFaces,
                expandedPoints
            );
        }
    }
    else
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcExpandedSlave() const")
            << "Attempting to create expanded slave on a shadow"
            << abort(FatalError);
    }
}


const Foam::standAlonePatch& Foam::overlapGgiPolyPatch::expandedMaster() const
{
    if (!expandedMasterPtr_)
    {
        calcExpandedMaster();
    }

    return *expandedMasterPtr_;
}

const Foam::standAlonePatch& Foam::overlapGgiPolyPatch::expandedSlave() const
{
    if (!expandedSlavePtr_)
    {
        calcExpandedSlave();
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
                0              // slave overlap tolerance
            );
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
    // Create neighbouring face centres using interpolation
    if (master())
    {
        const label shadowID = shadowIndex();

        // Get the transformed and interpolated shadow face cell centers
        vectorField delta = boundaryMesh()[shadowID].faceCellCentres()
            - boundaryMesh()[shadowID].faceCentres();

        reconFaceCellCentresPtr_ =
            new vectorField(interpolate(delta) + faceCentres());
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


void Foam::overlapGgiPolyPatch::clearOut()
{
    deleteDemandDrivenData(expandedMasterPtr_);
    deleteDemandDrivenData(expandedSlavePtr_);
    deleteDemandDrivenData(patchToPatchPtr_);
    deleteDemandDrivenData(reconFaceCellCentresPtr_);
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


void Foam::overlapGgiPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}


void Foam::overlapGgiPolyPatch::calcGeometry()
{
    // Reconstruct the cell face centres
    if (patchToPatchPtr_ && master())
    {
        reconFaceCellCentres();
    }

    calcTransforms();
    polyPatch::calcGeometry();
}


void Foam::overlapGgiPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}


void Foam::overlapGgiPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);

    // Force recalculation of interpolation
    clearOut();
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
