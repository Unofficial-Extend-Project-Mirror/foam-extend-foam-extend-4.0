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

\*---------------------------------------------------------------------------*/

#include "twoStrokeEngine.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoStrokeEngine::makeLayersLive()
{
    const polyTopoChanger& morphs = topoChanger_;

    // Enable layering
    forAll (morphs, modI)
    {
        if (isA<layerAdditionRemoval>(morphs[modI]))
        {
            morphs[modI].enable();
        }
        else if (isA<slidingInterface>(morphs[modI]))
        {
            morphs[modI].disable();
        }
        else
        {
            FatalErrorIn("void Foam::twoStrokeEngine::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << morphs[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::twoStrokeEngine::makeSlidersLive()
{
    const polyTopoChanger& morphs = topoChanger_;

    // Enable sliding interface
    forAll (morphs, modI)
    {
        if (isA<layerAdditionRemoval>(morphs[modI]))
        {
            morphs[modI].disable();
        }
        else if (isA<slidingInterface>(morphs[modI]))
        {
            morphs[modI].enable();
        }
        else
        {
            FatalErrorIn("void movingSquaresTM::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << morphs[modI].type()
                << abort(FatalError);
        }
    }

}

bool Foam::twoStrokeEngine::attached() const
{
    const polyTopoChanger& morphs = topoChanger_;

    bool result = false;

    forAll (morphs, modI)
    {
        if (isA<slidingInterface>(morphs[modI]))
        {
            result = result
             || refCast<const slidingInterface>(morphs[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll (morphs, modI)
    {
        if (isA<slidingInterface>(morphs[modI]))
        {
            if
            (
                result
             != refCast<const slidingInterface>(morphs[modI]).attached()
            )
            {
                FatalErrorIn("bool movingSquaresTM::attached() const")
                    << "Slider " << modI << " named " << morphs[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    // Sync across processors
    reduce(result, orOp<bool>());

    return result;
}


bool Foam::twoStrokeEngine::update()
{
    // Detaching the interface
    if (attached())
    {
        Info << "Decoupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

        bool localMorphing1 = topoChangeMap1->morphing();

        // Note: Since we are detaching, global morphing is always true
        // HJ, 7/Mar/2011

        if (localMorphing1)
        {
            Info << "Topology change; executing pre-motion after "
                << "sliding detach" << endl;
            movePoints(topoChangeMap1->preMotionPoints());
        }
        else
        {
            pointField newPoints = allPoints();

            // Dummy motion
            movePoints(newPoints);
        }

        Info << "sliding interfaces successfully decoupled!!!" << endl;
    }
    else
    {
        Info << "Sliding interfaces decoupled" << endl;
    }

    Info << "Executing layer action" << endl;

    // Piston Layering

    makeLayersLive();

    // Find piston mesh modifier if present on processor
    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    scalar minLayerThickness = piston().minLayer();
    scalar deltaZ = engTime().pistonDisplacement().value();
    virtualPistonPosition_ += deltaZ;

    Info << "virtualPistonPosition = " << virtualPistonPosition()
    << ", deckHeight = " << deckHeight()
    << ", pistonPosition = " << pistonPosition() << endl;

    if (realDeformation())
    {
        // Dectivate piston layer
        if (pistonLayerID > -1)
        {
            Info << "Mesh deformation mode" << endl;
            topoChanger_[pistonLayerID].disable();
        }
    }
    else
    {
        // Activate piston layer
        if (pistonLayerID > -1)
        {
            Info << "Piston layering mode" << endl;
            topoChanger_[pistonLayerID].enable();
        }
    }


    // Changing topology by hand
    autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

    bool localMorphing2 = topoChangeMap2->morphing();
    bool globalMorphing2 = localMorphing2;

    //HJ Is this reduce missing?  HJ, 19/Nov/2013
    reduce(globalMorphing2, orOp<bool>());

    // Work array for new points position.
    pointField newPoints = allPoints();

    if (globalMorphing2)
    {
        Info<< "Topology change; executing pre-motion after "
            << "dynamic layering" << endl;

        if (localMorphing2)
        {
            movePoints(topoChangeMap2->preMotionPoints());
            newPoints = topoChangeMap2->preMotionPoints();
        }
        else
        {
            movePoints(newPoints);
        }

        setV0();
        resetMotion();
    }

    // Reset the position of layered interfaces

    boolList scaleDisp(nPoints(), true);

    boolList pistonPoint(newPoints.size(), false);
    boolList headPoint(newPoints.size(), false);

    // Make a list of piston and head points. HJ, 2/Dec/2009

    labelList pistonPoints;
    labelList headPoints;
    {
        // Get cell-point addressing
        const labelListList& cp = cellPoints();

        boolList count(newPoints.size(), false);

        // Piston points
        label pistonCellID = cellZones().findZoneID("pistonCells");

        if (pistonCellID > -1)
        {
            const labelList& pistonCells = cellZones()[pistonCellID];

            forAll (pistonCells, cellI)
            {
                const labelList& curCellPoints = cp[pistonCells[cellI]];

                forAll (curCellPoints, i)
                {
                    count[curCellPoints[i]] = true;
                }
            }

            // Count the points
            label nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    nCounted++;
                }
            }

            pistonPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    pistonPoints[nCounted] = pointI;
                    nCounted++;
                }
            }
        }

        // Repeat for head points
        count = false;

        const label headCellID = cellZones().findZoneID("headCells");

        if (headCellID > -1)
        {
            const labelList& headCells = cellZones()[headCellID];

            forAll (headCells, cellI)
            {
                const labelList& curCellPoints = cp[headCells[cellI]];

                forAll (curCellPoints, i)
                {
                    count[curCellPoints[i]] = true;
                }
            }

            // Count the points
            label nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    nCounted++;
                }
            }

            headPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    headPoints[nCounted] = pointI;
                    nCounted++;
                }
            }
        }
    }


    label nScaled = nPoints();

    const scalarField& movingPointsM = movingPointsMask();

    forAll (pistonPoints, i)
    {
        label pointI = pistonPoints[i];
        pistonPoint[pointI] = true;
        point& p = newPoints[pointI];

        if (p.z() < pistonPosition() - 1.0e-6)
        {
            scaleDisp[pointI] = false;
            nScaled--;
        }
    }

    forAll (headPoints, i)
    {
        headPoint[headPoints[i]] = true;
        scaleDisp[headPoints[i]] = false;
        nScaled--;
    }

    if (realDeformation())
    {
        forAll (scaleDisp, pointI)
        {
            point& p = newPoints[pointI];

            if (scaleDisp[pointI])
            {
                p.z() += movingPointsM[pointI]*deltaZ*
                    (deckHeight() - p.z())/(deckHeight() - pistonPosition());
            }
            else
            {
                if (!headPoint[pointI])
                {
                    p.z() += movingPointsM[pointI]*deltaZ;
                }
            }
        }
    }
    else
    {
        // Always move piston
        scalar pistonTopZ = -GREAT;

        forAll (pistonPoints, i)
        {
            point& p = newPoints[pistonPoints[i]];
            p.z() += deltaZ*movingPointsM[pistonPoints[i]];
            pistonTopZ = max(pistonTopZ, p.z());
        }


        // NN! fix. only needed for compression
        if (deltaZ > 0.0)
        {
            // check if piston-points have moved beyond the layer above
            forAll (newPoints, i)
            {
                if (!pistonPoint[i])
                {
                    if (virtualPistonPosition() > newPoints[i].z())
                    {
                        newPoints[i].z() =
                        (1.0 - movingPointsM[i])*newPoints[i].z()
                        +
                        movingPointsM[i] *
                        (pistonTopZ + 0.9*minLayerThickness);
                    }
                }
            }
        }
    }

    movePoints(newPoints);
    deleteDemandDrivenData(movingPointsMaskPtr_);

    pistonPosition() += deltaZ;

    // Changing topology by hand
    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        // Attach the interface
        Info << "Coupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        bool localMorphing3 = topoChangeMap3->morphing();
        bool globalMorphing3 = localMorphing3;

        reduce(globalMorphing3, orOp<bool>());

        if (globalMorphing3)
        {
            Info<< "Topology change; executing pre-motion after "
                << "sliding attach" << endl;

            // Grab points
            newPoints = allPoints();

            if (localMorphing3)
            {
                // If there is layering, pick up correct points
                if (topoChangeMap3->hasMotionPoints())
                {
                    newPoints = topoChangeMap3->preMotionPoints();
                }

                // Prepare old points for the move
                pointField mappedOldPointsNew(allPoints().size());

                mappedOldPointsNew.map
                (
                    oldPointsNew, topoChangeMap3->pointMap()
                );

                forAll (scavInPortPatches_, patchi)
                {
                    // Find cut point zone ID
                    const label cutPointZoneID = pointZones().findZoneID
                    (
                        "cutPointZone" + Foam::name(patchi + 1)
                    );

                    if (cutPointZoneID > -1)
                    {
                        const labelList& cutPointsAddressing =
                            pointZones()[cutPointZoneID];

                        forAll (cutPointsAddressing, i)
                        {
                            mappedOldPointsNew[cutPointsAddressing[i]] =
                                newPoints[cutPointsAddressing[i]];
                        }

                        forAll (cutPointsAddressing, i)
                        {
                            if
                            (
                                newPoints[cutPointsAddressing[i]].z()
                              > virtualPistonPosition()
                            )
                            {
                                mappedOldPointsNew
                                    [cutPointsAddressing[i]].z() =
                                    newPoints[cutPointsAddressing[i]].z();
                            }
                            else
                            {
                                mappedOldPointsNew
                                    [cutPointsAddressing[i]].z() =
                                    newPoints[cutPointsAddressing[i]].z()
                                  - deltaZ;
                            }
                        }
                    }
                }

                // Move mesh into correct old configuration
                movePoints(mappedOldPointsNew);

                resetMotion();
                setV0();

                // Set new point motion
                movePoints(newPoints);
            }
            else
            {
                // No local topological change.  Execute double motion for
                // sync with topological changes
                movePoints(oldPointsNew);

                resetMotion();
                setV0();

                // Set new point motion
                movePoints(newPoints);
            }
        }
    }

    return true;
}


// ************************************************************************* //
