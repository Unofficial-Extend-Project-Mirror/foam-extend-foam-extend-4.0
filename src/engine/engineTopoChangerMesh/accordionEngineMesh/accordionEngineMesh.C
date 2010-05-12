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

#include "accordionEngineMesh.H"
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
#include "symmetryFvPatch.H"
#include "tetDecompositionMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(accordionEngineMesh, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, accordionEngineMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



    
bool Foam::accordionEngineMesh::realDeformation() const
{

    if (virtualPistonPosition() + engTime().pistonDisplacement().value() > deckHeight_ - SMALL)
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
Foam::accordionEngineMesh::accordionEngineMesh
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    valves_(*this, engTime().engineDict().lookup("accordionEngineMesh")),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
    delta_(readScalar(engTime().engineDict().lookup("delta"))),
    offSet_(readScalar(engTime().engineDict().lookup("offSet"))),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    deckHeight_(GREAT),
    msPtr_(motionSolver::New(*this)),
    cylinderHeadName_(engTime().engineDict().lookup("cylinderHeadName")),
    linerName_(engTime().engineDict().lookup("linerName")),
    pistonAuxPoints_(engTime().engineDict().lookup("pistonAuxPoints")),
    moveDetach_(engTime().engineDict().lookup("moveDetach"))

{
    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::accordionEngineMesh::setBoundaryVelocity(volVectorField& U)
{
    // Set valve velociaty
    forAll (valves(), valveI)
    {
        vector valveVel =
            valves()[valveI].curVelocity()*valves()[valveI].cs().axis();


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
            U.boundaryField()[valves()[valveI].detachInPortPatchID().index()] ==
                vector::zero;
            U.oldTime().boundaryField()[valves()[valveI].detachInPortPatchID().index()] ==
                vector::zero;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].detachInCylinderPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].detachInCylinderPatchID().index()] ==
               vector::zero;
            U.oldTime().boundaryField()[valves()[valveI].detachInCylinderPatchID().index()] ==
               vector::zero;
        }
     
    }

}




// ************************************************************************* //
