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
        if (typeid(morphs[modI]) == typeid(layerAdditionRemoval))
        {
            morphs[modI].enable();
        }
        else if (typeid(morphs[modI]) == typeid(slidingInterface))
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
        if (typeid(morphs[modI]) == typeid(layerAdditionRemoval))
        {
            morphs[modI].disable();
        }
        else if (typeid(morphs[modI]) == typeid(slidingInterface))
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
        if (typeid(morphs[modI]) == typeid(slidingInterface))
        {
            result =
                result
             || refCast<const slidingInterface>(morphs[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll (morphs, modI)
    {
        if (typeid(morphs[modI]) == typeid(slidingInterface))
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
        autoPtr<mapPolyMesh> topoChangeMap5 = topoChanger_.changeMesh();

        if (topoChangeMap5->hasMotionPoints() && topoChangeMap5->morphing())
        {
            Info << "Topology change; executing pre-motion after "
                << "sliding detach" << endl;
            movePoints(topoChangeMap5->preMotionPoints());
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
    // Changing topology by hand


// /* Tommaso, 23/5/2008

    // Find piston mesh modifier
    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    if (pistonLayerID < 0)
    {
        FatalErrorIn("void engineFvMesh::moveAndMorph()")
            << "Piston modifier not found."
            << abort(FatalError);
    }

    scalar minLayerThickness = piston().minLayer();
    scalar deltaZ = engTime().pistonDisplacement().value();
    virtualPistonPosition() += deltaZ;

    Info << "virtualPistonPosition = " << virtualPistonPosition()
    << ", deckHeight = " << deckHeight()
    << ", pistonPosition = " << pistonPosition() << endl;

    if (realDeformation())
    {
        // Dectivate piston layer
        Info << "Mesh deformation mode" << endl;
        topoChanger_[pistonLayerID].disable();
    }
    else
    {
        // Activate piston layer
        Info << "Piston layering mode" << endl;
        topoChanger_[pistonLayerID].enable();
    }


    // Changing topology by hand
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    // Work array for new points position.
    pointField newPoints = allPoints();
    const pointField& refPoints = allPoints();

    if (topoChangeMap->morphing())
    {
        if (topoChangeMap->hasMotionPoints())
        {
            Info<< "Topology change; executing pre-motion after "
                << "dynamic layering" << endl;
            movePoints(topoChangeMap->preMotionPoints());
            newPoints = topoChangeMap->preMotionPoints();
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
        label pistonCellIndex = cellZones().findZoneID("pistonCells");

        if (pistonCellIndex < 0)
        {
            FatalErrorIn("bool twoStrokeEngine::update()")
                << "Cannot find cell zone pistonCells"
                << abort(FatalError);
        }


        const labelList& pistonCells = cellZones()[pistonCellIndex];

        label headCellIndex = cellZones().findZoneID("headCells");

        if (headCellIndex < 0)
        {
            FatalErrorIn("bool twoStrokeEngine::update()")
                << "Cannot find cell zone headCells"
                << abort(FatalError);
        }

        const labelList& headCells = cellZones()[headCellIndex];

        const labelListList& cp = cellPoints();

        boolList count(newPoints.size(), false);

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

        // Repeat for head points
        count = false;

        forAll (headCells, cellI)
        {
            const labelList& curCellPoints = cp[pistonCells[cellI]];

            forAll (curCellPoints, i)
            {
                count[curCellPoints[i]] = true;
            }
        }

        // Count the points
        nCounted = 0;
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


    label nScaled = nPoints();

//     label pistonPtsIndex = pointZones().findZoneID("pistonPoints");
//     const labelList& pistonPoints = pointZones()[pistonPtsIndex];

//     label headPtsIndex = pointZones().findZoneID("headPoints");
//     const labelList& headPoints = pointZones()[headPtsIndex];

    const scalarField& movingPointsM = movingPointsMask();

    forAll(pistonPoints, i)
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

    forAll(headPoints, i)
    {
        headPoint[headPoints[i]] = true;
        scaleDisp[headPoints[i]] = false;
        nScaled--;
    }

    if (realDeformation())
    {
        forAll(scaleDisp, pointI)
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
        forAll(pistonPoints, i)
        {
            point& p = newPoints[pistonPoints[i]];
            p.z() += deltaZ*movingPointsM[pistonPoints[i]];
            pistonTopZ = max(pistonTopZ, p.z());
        }


        // NN! fix. only needed for compression
        if (deltaZ > 0.0)
        {
            // check if piston-points have moved beyond the layer above
            forAll(newPoints, i)
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

//*/ //Tommaso, 23/5/2008

    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        // Attach the interface
        Info << "Coupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap4 = topoChanger_.changeMesh();

        if (topoChangeMap4->morphing())
        {
            if (topoChangeMap4->hasMotionPoints())
            {
                Info<< "Topology change; executing pre-motion after "
                    << "sliding attach" << endl;

//                 Info<< "topoChangeMap4->preMotionPoints().size() = "
//                     << topoChangeMap4->preMotionPoints().size() << nl
//                     << "allPoints.size() = " << allPoints().size() << nl
//                     << "points.size() = " << points().size() << endl;

                movePoints(topoChangeMap4->preMotionPoints());
                newPoints = points();
            }

            {
                pointField mappedOldPointsNew(allPoints().size());

                mappedOldPointsNew.map
                (
                    oldPointsNew, topoChangeMap4->pointMap()
                );

                forAll(scavInPortPatches_, patchi)
                {
                    const labelList& cutPointsAddressing =
                        pointZones()
                        [
                            pointZones().findZoneID
                            (
                                "cutPointZone" + Foam::name(patchi + 1)
                            )
                        ];

                    forAll(cutPointsAddressing, i)
                    {
                        mappedOldPointsNew[cutPointsAddressing[i]] =
                            newPoints[cutPointsAddressing[i]];
                    }

                    forAll(cutPointsAddressing, i)
                    {
                        if
                        (
                            newPoints[cutPointsAddressing[i]].z()
                          > virtualPistonPosition()
                        )
                        {
                            mappedOldPointsNew[cutPointsAddressing[i]].z() =
                                newPoints[cutPointsAddressing[i]].z();
                        }
                        else
                        {
                            mappedOldPointsNew[cutPointsAddressing[i]].z() =
                                newPoints[cutPointsAddressing[i]].z() - deltaZ;
                        }
                    }
                }
                pointField newPoints = allPoints();

                movePoints(mappedOldPointsNew);

                resetMotion();
                setV0();
                movePoints(newPoints);
            }
        }
    }

    return true;
}


// ************************************************************************* //
