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

#include "polyTopoChange.H"
#include "polyMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::polyTopoChange::minListSize = 100;
const Foam::label Foam::polyTopoChange::pointFraction = 10;
const Foam::label Foam::polyTopoChange::faceFraction = 10;
const Foam::label Foam::polyTopoChange::cellFraction = 10;

int Foam::polyTopoChange::debug
(
    Foam::debug::debugSwitch("polyTopoChange", 0)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh reference
Foam::polyTopoChange::polyTopoChange(const polyMesh& mesh)
:
    mesh_(mesh),
    addedPoints_(minListSize),
    modifiedPoints_(minListSize),
    removedPoints_
    (
        Foam::max(mesh.allPoints().size()/pointFraction, minListSize)
    ),
    addedFaces_(minListSize),
    modifiedFaces_(minListSize),
    removedFaces_(Foam::max(mesh.allFaces().size()/faceFraction, minListSize)),
    addedCells_(minListSize),
    removedCells_(Foam::max(mesh.nCells()/cellFraction, minListSize))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyTopoChange::~polyTopoChange()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::polyTopoChange::setAction(const topoAction& action)
{
    if (isType<polyAddPoint>(action))
    {
        addedPoints_.append(refCast<const polyAddPoint>(action));

        return mesh_.allPoints().size() + addedPoints_.size() - 1;
    }
    else if (isType<polyModifyPoint>(action))
    {
        const polyModifyPoint& pmp = refCast<const polyModifyPoint>(action);

        if (debug)
        {
            if (removedPoints_.find(pmp.pointID()) != removedPoints_.end())
            {
                FatalErrorIn
                (
                    "label polyTopoChange::setAction(const topoAction& action)"
                )   << "Modifying point " << pmp.pointID()
                    << " for the second time or modifying a removed point.  "
                    << "This is not allowed."
                    << abort(FatalError);
            }
        }

        modifiedPoints_.append(pmp);
        removedPoints_.insert(pmp.pointID());

        return -1;
    }
    else if (isType<polyRemovePoint>(action))
    {
        const polyRemovePoint& prp = refCast<const polyRemovePoint>(action);

        if (debug)
        {
            if (removedPoints_.find(prp.pointID()) != removedPoints_.end())
            {
                FatalErrorIn
                (
                    "label polyTopoChange::setAction(const topoAction& action)"
                )   << "Removing point " << prp.pointID()
                    << " for the second time or removing a modified point.  "
                    << "This is not allowed."
                    << abort(FatalError);
            }
        }

        removedPoints_.insert(prp.pointID());

        return -1;
    }
    else if (isType<polyAddFace>(action))
    {
        addedFaces_.append(refCast<const polyAddFace>(action));

        return mesh_.allFaces().size() + addedFaces_.size() - 1;
    }
    else if (isType<polyModifyFace>(action))
    {
        const polyModifyFace& pmf = refCast<const polyModifyFace>(action);

        if (debug)
        {
            if (removedFaces_.find(pmf.faceID()) != removedFaces_.end())
            {
                FatalErrorIn
                (
                    "label polyTopoChange::setAction(const topoAction& action)"
                )   << "Modifying face " << pmf.faceID()
                    << " for the second time or modifying a removed face.  "
                    << "This is not allowed."
                    << abort(FatalError);
            }
        }

        modifiedFaces_.append(pmf);
        removedFaces_.insert(pmf.faceID());

        return -1;
    }
    else if (isType<polyRemoveFace>(action))
    {
        const polyRemoveFace& prf = refCast<const polyRemoveFace>(action);

        if (debug)
        {
            if (removedFaces_.find(prf.faceID()) != removedFaces_.end())
            {
                FatalErrorIn
                (
                    "label polyTopoChange::setAction(const topoAction& action)"
                )   << "Removing face " << prf.faceID()
                    << " for the second time or removing a modified face.  "
                    << "This is not allowed."
                    << abort(FatalError);
            }
        }

        removedFaces_.insert(prf.faceID());

        return -1;
    }
    else if (isType<polyAddCell>(action))
    {
        addedCells_.append(refCast<const polyAddCell>(action));

        return mesh_.nCells() + addedCells_.size() - 1;
    }
    else if (isType<polyModifyCell>(action))
    {
        const polyModifyCell& pmc = refCast<const polyModifyCell>(action);

        modifiedCells_.append(pmc);

        return -1;
    }
    else if (isType<polyRemoveCell>(action))
    {
        const polyRemoveCell& prc = refCast<const polyRemoveCell>(action);

        if (debug)
        {
            if (removedCells_.find(prc.cellID()) != removedCells_.end())
            {
                FatalErrorIn
                (
                    "label polyTopoChange::setAction(const topoAction& action)"
                )   << "Removing cell " << prc.cellID()
                    << " for the second time.  This is not allowed."
                    << abort(FatalError);
            }
        }

        removedCells_.insert(prc.cellID());

        return -1;
    }
    else
    {
        FatalErrorIn
        (
            "label polyTopoChange::setAction(const topoAction& action)"
        )   << "Unknown type of topoChange: " << action.type()
            << abort(FatalError);

        // Dummy return to keep compiler happy
        return -1;
    }
}


bool Foam::polyTopoChange::check() const
{
    if (debug)
    {
        Info<< "Points: added = " << addedPoints_.size()
            << " modified = " << modifiedPoints_.size()
            << " removed = " << removedPoints_.size() - modifiedPoints_.size()
            << " balance = " << pointBalance()
            << endl;

        Info<< "Faces: added = " << addedFaces_.size()
            << " modified = " << modifiedFaces_.size()
            << " removed = " << removedFaces_.size() - modifiedFaces_.size()
            << " balance = " << faceBalance()
            << endl;

        Info<< "Cells: added = " << addedCells_.size()
            << " modified = " << modifiedCells_.size()
            << " removed = " << removedCells_.size() - modifiedCells_.size()
            << " balance = " << cellBalance()
            << endl;
    }

    label nErrors = 0;

    // Check points

    if (mesh_.allPoints().size() + pointBalance() < 4)
    {
        WarningIn("bool polyTopoChange::check() const")
            << "Less than 4 points in the mesh.  This cannot be a valid mesh."
            << endl;

        nErrors++;
    }

    // Collect zone-only and removed points
    labelHashSet zoneOnlyPoints(removedPoints_);

    forAll (modifiedPoints_, mpI)
    {
        if (modifiedPoints_[mpI].inCell())
        {
            // Point is live; remove from lookup
            zoneOnlyPoints.erase
            (
                zoneOnlyPoints.find(modifiedPoints_[mpI].pointID())
            );
        }
    }

    forAll (addedPoints_, apI)
    {
        if (!addedPoints_[apI].inCell())
        {
            zoneOnlyPoints.insert(mesh_.allPoints().size() + apI);
        }
    }

    // Check faces

    label nFaceErrors = 0;
    label nCurFaceError = 0;

    // For all live modified and added faces, check that points are also live

    // Collect zone-only faces
    labelHashSet zoneOnlyFaces(removedFaces_);

    // Modified faces
    forAll (modifiedFaces_, mfI)
    {
        if (!modifiedFaces_[mfI].onlyInZone())
        {
            // Face is live; remove from lookup
            zoneOnlyFaces.erase
            (
                zoneOnlyFaces.find(modifiedFaces_[mfI].faceID())
            );

            // Count errors in current face
            nCurFaceError = 0;

            const face& curFace = modifiedFaces_[mfI].newFace();

            forAll (curFace, pointI)
            {
                if (zoneOnlyPoints.found(curFace[pointI]))
                {
                    nCurFaceError++;
                }
            }

            if (nCurFaceError > 0)
            {
                if (debug)
                {
                    Info<< "Modified face no " << mfI << " is active but uses "
                        << nCurFaceError << " zone-only or removed points.  "
                        << "Face: " << curFace << "  Problem points:";

                    forAll (curFace, pointI)
                    {
                        if (zoneOnlyPoints.found(curFace[pointI]))
                        {
                            Info  << " " << curFace[pointI];
                        }
                    }

                    Info << endl;
                }

                nFaceErrors++;
            }
        }
    }

    // Added faces
    forAll (addedFaces_, afI)
    {
        if (addedFaces_[afI].onlyInZone())
        {
            // Zone-only face; add to hash set
            zoneOnlyFaces.insert(mesh_.allFaces().size() + afI);
        }
        else
        {
            // Count errors in current face
            nCurFaceError = 0;

            const face& curFace = addedFaces_[afI].newFace();

            forAll (curFace, pointI)
            {
                if (zoneOnlyPoints.found(curFace[pointI]))
                {
                    nCurFaceError++;
                }
            }

            if (nCurFaceError > 0)
            {
                if (debug)
                {
                    Info<< "Added face no " << afI << " is active but uses "
                        << nCurFaceError << " zone-only or removed points.  "
                        << "Face: " << curFace << "  Problem points:";

                    forAll (curFace, pointI)
                    {
                        if (zoneOnlyPoints.found(curFace[pointI]))
                        {
                            Info  << " " << curFace[pointI];
                        }
                    }

                    Info << endl;
                }

                nFaceErrors++;
            }
        }
    }

    // Check cells?  

    if (nErrors == 0)
    {
        if (debug)
        {
            Info << "polyTopoChange OK" << endl;
        }

        return false;
    }
    else
    {
        if (debug)
        {
            Info<< "bool polyTopoChange::check() const : "
                << "Detected " << nErrors << " errors."
                << endl;
        }

        return true;
    }
}


void Foam::polyTopoChange::clear()
{
    if (debug)
    {
        Info<< "void Foam::polyTopoChange::clear() : "
            << "clearing topological change request" << endl;
    }

    // Clear all contents
    addedPoints_.clear();
    modifiedPoints_.clear();
    removedPoints_.clear();

    addedFaces_.clear();
    modifiedFaces_.clear();
    removedFaces_.clear();

    addedCells_.clear();
    removedCells_.clear();
}


// ************************************************************************* //
