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

\*---------------------------------------------------------------------------*/

#include "RBFMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBFMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        RBFMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBFMotionSolver::makeControlIDs()
{
    // Points that are neither on moving nor on static patches
    // will be marked with 0
    labelList markedPoints(mesh().nPoints(), 0);

    // Mark all points on moving patches with 1
    label nMarkedPoints = 0;

    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(movingPatches_[patchI]);

        if (patchIndex < 0)
        {
            FatalErrorIn("void RBFMotionSolver::makeControlIDs()")
                << "Patch " << movingPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh().boundaryMesh().names()
                << abort(FatalError);
        }

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
            markedPoints[mp[i]] = 1;
            nMarkedPoints++;
        }
    }

    // Mark moving points and select control points from moving patches
    movingIDs_.setSize(nMarkedPoints);
    controlIDs_.setSize(nMarkedPoints);

    Info << "Total points on moving boundaries: " << nMarkedPoints << endl;

    const pointField& points = mesh().points();

    // Re-use counter
    nMarkedPoints = 0;

    forAll (markedPoints, i)
    {
        if (markedPoints[i] == 1)
        {
            // Grab internal point
            movingIDs_[nMarkedPoints] = i;
            nMarkedPoints++;
        }
    }

    movingIDs_.setSize(nMarkedPoints);

    // Actual location of moving points will be set later on request
    // HJ, 19/Dec/2008
    movingPoints_.setSize(nMarkedPoints, vector::zero);
    motion_.setSize(nMarkedPoints, vector::zero);

    // Re-use counter
    nMarkedPoints = 0;

    // Mark all points on static patches with -1
    forAll (staticPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(staticPatches_[patchI]);

        if (patchIndex < 0)
        {
            FatalErrorIn("void RBFMotionSolver::makeControlPoints()")
                << "Patch " << staticPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh().boundaryMesh().names()
                << abort(FatalError);
        }

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
            markedPoints[mp[i]] = -1;
            nMarkedPoints++;
        }
    }

    Info << "Total points on static boundaries: " << nMarkedPoints << endl;

    // Re-use counter
    nMarkedPoints = 0;

    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(movingPatches_[patchI]);

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        for
        (
            label pickedPoint = 0;
            pickedPoint < mp.size();
            pickedPoint += coarseningRatio_
        )
        {
            // Pick point as control point
            controlIDs_[nMarkedPoints] = mp[pickedPoint];

            // Mark the point as picked
            markedPoints[mp[pickedPoint]] = 2;
            nMarkedPoints++;
        }
    }

    Info << "Selected " << nMarkedPoints << " control points" << endl;

    // Resize control IDs
    controlIDs_.setSize(nMarkedPoints);

    // Pick up point locations
    controlPoints_.setSize(nMarkedPoints);

    // Set control points
    forAll (controlIDs_, i)
    {
        controlPoints_[i] = points[controlIDs_[i]];
    }

    // Pick up all internal points
    internalIDs_.setSize(points.size());
    internalPoints_.setSize(points.size());

    // Re-use counter
    nMarkedPoints = 0;

    forAll (markedPoints, i)
    {
        if (markedPoints[i] == 0)
        {
            // Grab internal point
            internalIDs_[nMarkedPoints] = i;
            internalPoints_[nMarkedPoints] = points[i];
            nMarkedPoints++;
        }
    }

    Info << "Number of internal points: " << nMarkedPoints << endl;

    // Resize the lists
    internalIDs_.setSize(nMarkedPoints);
    internalPoints_.setSize(nMarkedPoints);
}


void Foam::RBFMotionSolver::setMovingPoints() const
{
    const pointField& points = mesh().points();

    // Set moving points
    forAll (movingIDs_, i)
    {
        movingPoints_[i] = points[movingIDs_[i]];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFMotionSolver::RBFMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    motionSolver(mesh),
    movingPatches_(lookup("movingPatches")),
    staticPatches_(lookup("staticPatches")),
    coarseningRatio_(readLabel(lookup("coarseningRatio"))),
    movingIDs_(0),
    movingPoints_(0),
    controlIDs_(0),
    controlPoints_(0),
    internalIDs_(0),
    internalPoints_(0),
    motion_(0),
    interpolation_
    (
        subDict("interpolation"),
        controlPoints_,
        internalPoints_
    )
{
    makeControlIDs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBFMotionSolver::~RBFMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBFMotionSolver::setMotion(const vectorField& m)
{
    motion_ = m;

    // Re-calculate interpolation
    interpolation_.movePoints();
}


const Foam::vectorField& Foam::RBFMotionSolver::movingPoints() const
{
    // Update moving points based on current mesh
    setMovingPoints();

    return movingPoints_;
}


Foam::tmp<Foam::pointField> Foam::RBFMotionSolver::curPoints() const
{
    // Prepare new points: same as old point
    tmp<pointField> tnewPoints
    (
        new vectorField(mesh().nPoints(), vector::zero)
    );
    pointField& newPoints = tnewPoints();

    // Add motion to existing points

    // 1. Insert prescribed motion of moving points
    forAll (movingIDs_, i)
    {
        newPoints[movingIDs_[i]] = motion_[i];
    }

    vectorField motionOfControl(controlIDs_.size());

    // 2. Capture positions of control points
    forAll (controlIDs_, i)
    {
        motionOfControl[i] = newPoints[controlIDs_[i]];
    }

    // Call interpolation
    vectorField interpolatedMotion =
        interpolation_.interpolate(motionOfControl);

    // 3. Insert RBF interpolated motion
    forAll (internalIDs_, i)
    {
        newPoints[internalIDs_[i]] = interpolatedMotion[i];
    }

    // 4. Add old point positions
    newPoints += mesh().points();

    return tnewPoints;
}


void Foam::RBFMotionSolver::solve()
{}


void Foam::RBFMotionSolver::updateMesh(const mapPolyMesh&)
{
    // Recalculate control point IDs
    makeControlIDs();
}


// ************************************************************************* //
