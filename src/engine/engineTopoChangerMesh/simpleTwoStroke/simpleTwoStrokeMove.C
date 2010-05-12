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


#include "simpleTwoStroke.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::simpleTwoStroke::makeLayersLive()
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
            FatalErrorIn("void Foam::simpleTwoStroke::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << morphs[modI].type()
                << abort(FatalError);
        }
    }
}

void Foam::simpleTwoStroke::makeSlidersLive()
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

bool Foam::simpleTwoStroke::attached() const
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


bool Foam::simpleTwoStroke::update()
{
    // Detaching the interface
    if (attached())
    {
        Info << "Decoupling sliding interfaces" << endl;
        makeSlidersLive();
        topoChanger_.changeMesh();

        Info << "sliding interfaces successfully decoupled!!!" << endl;
    }
    else
    {
        Info << "Sliding interfaces decoupled" << endl;
    }

    Info << "Executing layer action" << endl;

    // Piston Layering

    makeLayersLive();

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
    pointField newPoints = points();

    if (topoChangeMap->morphing())
    {

        if (topoChangeMap->hasMotionPoints())
        {
            Info << "Topology change; executing pre-motion" << endl;
            movePoints(topoChangeMap->preMotionPoints());
            newPoints = topoChangeMap->preMotionPoints();
        }
        setV0();
        resetMotion();
    }

    // Reset the position of layered interfaces

    boolList scaleDisp(nPoints(), true);
    label nScaled = nPoints();

    List<bool> pistonPoint(newPoints.size(), false);
    List<bool> headPoint(newPoints.size(), false);

//    label pistonPtsIndex = pointZones().findZoneID("pistonPoints");
//    const labelList& pistonPoints = pointZones()[pistonPtsIndex];
    
    labelList pistonPoints;

    {
        label movingCellsIndex = cellZones().findZoneID("movingCells");

        if (movingCellsIndex < 0)
        {
            FatalErrorIn("bool twoStrokeEngine::update()")
                << "Cannot find cell zone movingCells"
                << abort(FatalError);
        }


        const labelList& pistonCells = cellZones()[movingCellsIndex];

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

    }

    label headPtsIndex = pointZones().findZoneID("headPoints");
    const labelList& headPoints = pointZones()[headPtsIndex];

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
                p.z() += movingPointsM[pointI]*
                    deltaZ
                  * (deckHeight() - p.z())/(deckHeight() - pistonPosition());
            }
            else
            {
                if (!headPoint[pointI])
                {
                    p.z() += movingPointsM[pointI] * deltaZ;
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
                          + movingPointsM[i]*
                            (pistonTopZ + 0.9*minLayerThickness);
                    }
                }
            }
        }
    }


    movePoints(newPoints);

    deleteDemandDrivenData(movingPointsMaskPtr_);

    pistonPosition() += deltaZ;

    {
        // Attach the interface
        Info << "Coupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        Info << "Sliding interfaces coupled: " << attached() << endl;
    }

    return true;
}

