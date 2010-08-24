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

#include "layerSmooth.H"
#include "layerAdditionRemoval.H"
#include "attachDetach.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "engineTime.H"
#include "pointField.H"
#include "fvPatchField.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::layerSmooth, 0);

addToRunTimeSelectionTable
(
    Foam::engineTopoChangerMesh,
    Foam::layerSmooth,
    IOobject
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::layerSmooth::realDeformation() const
{

    if
    (
        virtualPistonPosition() + engTime().pistonDisplacement().value()
      > deckHeight_ - SMALL
    )
    {
        return true;
    }
    else
    {
        return deformation();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::layerSmooth::layerSmooth
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    valves_(*this, engTime().engineDict().lookup("layerSmooth")),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
    delta_(readScalar(engTime().engineDict().lookup("delta"))),
    offSet_(readScalar(engTime().engineDict().lookup("offSet"))),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    deckHeight_(GREAT),
    msPtr_(motionSolver::New(*this)),
    valvePlaneTol_(readScalar(engTime().engineDict().lookup("valvePlaneTol")))
{


    forAll(valves(), valveI)
    {
        if
        (
            valves()[valveI].stemPatchID().active()
         && valves()[valveI].bottomPatchID().active()
        )
        {
            label bottomIndex = valves()[valveI].bottomPatchID().index();
            label stemIndex = valves()[valveI].stemPatchID().index();

            const polyPatch& stemPatch =
                    boundaryMesh()[stemIndex];

            const polyPatch& bottomPatch =
                    boundaryMesh()[bottomIndex];

            point bottomPatchCentre = vector::zero;
            point stemPatchCentre = vector::zero;

            scalar bottomArea = 0;
            scalar stemArea = 0;

            Switch halfGeometry(engTime().engineDict().lookup("halfGeometry"));

            forAll(stemPatch, i)
            {
                stemPatchCentre +=
                    stemPatch.faceCentres()[i]*mag(stemPatch.faceAreas()[i]);

                stemArea += mag(stemPatch.faceAreas()[i]);
            }

            forAll(bottomPatch, i)
            {
                bottomPatchCentre +=bottomPatch.faceCentres()[i]*
                    mag(bottomPatch.faceAreas()[i]);

                bottomArea += mag(bottomPatch.faceAreas()[i]);
            }

            if (halfGeometry)
            {
                Info << "half Geometry active" << endl;

                forAll(stemPatch, i)
                {
                    stemPatchCentre +=
                        vector
                        (
                            stemPatch.faceCentres()[i].x(),
                            -1.0*stemPatch.faceCentres()[i].y(),
                            stemPatch.faceCentres()[i].z()
                        )*mag(stemPatch.faceAreas()[i]);

                    stemArea += mag(stemPatch.faceAreas()[i]);
                }

                forAll(bottomPatch, i)
                {
                    bottomPatchCentre +=
                        vector
                        (
                            bottomPatch.faceCentres()[i].x(),
                            -1.0*bottomPatch.faceCentres()[i].y(),
                            bottomPatch.faceCentres()[i].z()
                        )*mag(bottomPatch.faceAreas()[i]);

                    bottomArea += mag(bottomPatch.faceAreas()[i]);
                }
            }

            stemPatchCentre /= stemArea;
            bottomPatchCentre /= bottomArea;


            Info<< "Foam::layerSmooth::layerSmooth"
                << " calculated origin for valve n. " << valveI
                << " is " <<  (stemPatchCentre) << endl;

            Info<< "Foam::layerSmooth::layerSmooth"
                << " calculated axis for valve n. " << valveI
                << " is "
                << (stemPatchCentre - bottomPatchCentre)/
                mag(stemPatchCentre-bottomPatchCentre)
                << endl;
        }
    }

    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::layerSmooth::setBoundaryVelocity(volVectorField& U)
{
    // Set valve velociaty
    forAll (valves(), valveI)
    {
        vector valveVel =
            valves()[valveI].curVelocity()*valves()[valveI].cs().axis();

        // If valve is present in geometry, set the motion
        if (valves()[valveI].curtainInPortPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()
                [valves()[valveI].curtainInPortPatchID().index()] == valveVel;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].curtainInCylinderPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()
                [valves()[valveI].curtainInCylinderPatchID().index()] ==
                valveVel;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].stemPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].stemPatchID().index()] ==
                valveVel;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].detachInPortPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()
                [valves()[valveI].detachInPortPatchID().index()] ==
                vector::zero;

            U.oldTime().boundaryField()
                [valves()[valveI].detachInPortPatchID().index()] ==
                vector::zero;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].detachInCylinderPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()
                [valves()[valveI].detachInCylinderPatchID().index()] ==
               vector::zero;

            U.oldTime().boundaryField()
                [valves()[valveI].detachInCylinderPatchID().index()] ==
               vector::zero;
        }
    }
}


bool Foam::layerSmooth::isACylinderHeadFace
(
    const labelList& cylHeadFaces,
    const label face
)
{
    if(face >= cylHeadFaces[0] && face <= cylHeadFaces[cylHeadFaces.size()-1])
    {
        return true;
    }

/*
    forAll(cylHeadFaces, i)
    {
        if(cylHeadFaces[i] == face)
        {
            return true;
        }
    }
*/

    return false;
}



// ************************************************************************* //
