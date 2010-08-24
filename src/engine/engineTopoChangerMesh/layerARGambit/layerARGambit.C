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

#include "layerARGambit.H"
#include "engineTime.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"
#include "volPointInterpolation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(layerARGambit, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, layerARGambit, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layerARGambit::makeLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll (topoChanges, modI)
    {
        if (typeid(topoChanges[modI]) == typeid(layerAdditionRemoval))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorIn("void Foam::layerARGambit::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::layerARGambit::checkAndCalculate()
{
    label pistonIndex = -1;
    bool foundPiston = false;

    label linerIndex = -1;
    bool foundLiner = false;

    label cylinderHeadIndex = -1;
    bool foundCylinderHead = false;

    forAll(boundary(), i)
    {
        Info << boundary()[i].name() << endl;
        if (boundary()[i].name() == "piston")
        {
            pistonIndex = i;
            foundPiston = true;
        }
        else if (boundary()[i].name() == "liner")
        {
            linerIndex = i;
            foundLiner = true;
        }
        else if (boundary()[i].name() == "cylinderHead")
        {
            cylinderHeadIndex = i;
            foundCylinderHead = true;
        }
    }

    reduce(foundPiston, orOp<bool>());
    reduce(foundLiner, orOp<bool>());
    reduce(foundCylinderHead, orOp<bool>());

    if (!foundPiston)
    {
        FatalErrorIn("Foam::layerARGambit::checkAndCalculate()")
            << " : cannot find piston patch"
            << abort(FatalError);
    }

    if (!foundLiner)
    {
        FatalErrorIn("Foam::layerARGambit::checkAndCalculate()")
            << " : cannot find liner patch"
            << abort(FatalError);
    }

    if (!foundCylinderHead)
    {
        FatalErrorIn("Foam::layerARGambit::checkAndCalculate()")
            << " : cannot find cylinderHead patch"
            << exit(FatalError);
    }

    {
        if (linerIndex != -1)
        {
            pistonPosition() =
                max(boundary()[pistonIndex].patch().localPoints()).z();
        }
        reduce(pistonPosition(), minOp<scalar>());

        if (cylinderHeadIndex != -1)
        {
            deckHeight() = min
            (
                boundary()[cylinderHeadIndex].patch().localPoints()
            ).z();
        }
        reduce(deckHeight(), minOp<scalar>());

        Info<< "deckHeight: " << deckHeight() << nl
            << "piston position: " << pistonPosition() << endl;
    }
}

void Foam::layerARGambit::setVirtualPistonPosition()
{

    label pistonFaceIndex = faceZones().findZoneID("pistonLayerFaces");

    bool foundPistonFace = (pistonFaceIndex != -1);

    Info << "piston face index = " << pistonFaceIndex << endl;

    if(!foundPistonFace)
    {
        FatalErrorIn("Foam::layerARGambit::setVirtualPistonPosition()")
            << " : cannot find the pistonLayerFaces"
            << exit(FatalError);
    }

    const labelList& pistonFaces = faceZones()[pistonFaceIndex];
    forAll(pistonFaces, i)
    {
        const face& f = faces()[pistonFaces[i]];

        // should loop over facepoints...
        forAll(f, j)
        {
            virtualPistonPosition() =
                Foam::max(virtualPistonPosition(), points()[f[j]].z());
        }
    }

    reduce(virtualPistonPosition(), maxOp<scalar>());

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::layerARGambit::layerARGambit
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
    delta_(readScalar(engTime().engineDict().lookup("delta"))),
    offSet_(readScalar(engTime().engineDict().lookup("offSet"))),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    deckHeight_(GREAT)
{
    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::layerARGambit::~layerARGambit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::layerARGambit::update()
{
    scalar deltaZ = engTime().pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << " Piston at:" << pistonPosition()
        << endl;
    virtualPistonPosition() += deltaZ;

    // Find piston mesh modifier
    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    Info << "pistonLayerID: " << pistonLayerID << endl;

    const layerAdditionRemoval& pistonLayers =
        dynamicCast<const layerAdditionRemoval>
        (
            topoChanger_[pistonLayerID]
        );

    bool realDeformation = deformation();

    if
    (
        virtualPistonPosition() + deltaZ
     > deckHeight()-engTime().clearance().value() - SMALL
    )
    {
        realDeformation = true;
    }

    if (realDeformation)
    {
        // Disable layer addition
        Info << "Mesh deformation mode" << endl;
        topoChanger_[pistonLayerID].disable();
    }
    else
    {
        // Enable layer addition
        Info << "Piston layering mode" << endl;
        topoChanger_[pistonLayerID].enable();
    }

    scalar minLayerThickness = pistonLayers.minLayerThickness();

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    bool meshChanged = topoChangeMap->morphing();

    pointField newPoints = allPoints();

    if (meshChanged)
    {
        if (topoChangeMap().hasMotionPoints())
        {
            movePoints(topoChangeMap().preMotionPoints());
            newPoints = topoChangeMap().preMotionPoints();
        }
        setV0();
        resetMotion();
    }

    Info << "virtualPistonPosition = " << virtualPistonPosition()
    << ", deckHeight = " << deckHeight() << endl;

    // Mesh in three parts:
    // - pistonPoints - move with deltaZ
    // - headPoints - do not move

    const pointZoneMesh& pZones = pointZones();
    label headPtsIndex = pZones.findZoneID("headPoints");
    label pistonPtsIndex = pZones.findZoneID("pistonPoints");
    const labelList& pistonPoints = pZones[pistonPtsIndex];
    const labelList& headPoints = pZones[headPtsIndex];


    // Whether point displacement is by scaling
    boolList scaleDisp(nPoints(), true);
    label nScaled = nPoints();
    List<bool> pistonPoint(newPoints.size(), false);
    List<bool> headPoint(newPoints.size(), false);

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

    if (realDeformation)
    {

        forAll(scaleDisp, pointI)
        {
            point& p = newPoints[pointI];

            if (scaleDisp[pointI])
            {
                p.z() +=
                    deltaZ
                  * (deckHeight() - p.z())/(deckHeight() - pistonPosition());
            }
            else
            {
                if (!headPoint[pointI])
                {
                    p.z() += deltaZ;
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
            p.z() += deltaZ;
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
                        newPoints[i].z() = pistonTopZ + 0.9*minLayerThickness;
                    }
                }
            }
        }
    }

    movePoints(newPoints);

    pistonPosition() += deltaZ;
    scalar pistonSpeed = deltaZ/engTime().deltaT().value();

    Info<< "clearance: " << deckHeight() - pistonPosition() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;

    Info << "Total cylinder volume at CA " << engTime().timeName() << " = " <<
    sum(V()) << endl;

    return meshChanged;

}

void Foam::layerARGambit::setBoundaryVelocity(volVectorField& U)
{
// Does nothing, using the movingWallVelocity boundary condition for U in the piston patch...

//    vector pistonVel = piston().cs().axis()*engTime().pistonSpeed().value();

    //mean piston velocityy
/*
    vector pistonVel = 0.5 * piston().cs().axis()*
                            dimensionedScalar
                            (
                                "meanPistonSpeed",
                                dimLength/dimTime,
                                2.0*engTime().stroke().value()*engTime().rpm().value()/60.0
                            ).value();
*/

//    U.boundaryField()[piston().patchID().index()] = pistonVel;
}

// ************************************************************************* //
