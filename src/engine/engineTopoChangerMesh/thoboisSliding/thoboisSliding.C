/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2010 Tommaso Lucchini
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

Class
    verticalValvesGambit

\*---------------------------------------------------------------------------*/

#include "thoboisSliding.H"
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
    defineTypeNameAndDebug(thoboisSliding, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, thoboisSliding, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



    
bool Foam::thoboisSliding::realDeformation() const
{

    bool deformationValve = false;
    forAll(valves(), valveI)
    {
        
        scalar maxLayer = piston().minLayer();
            
        if(valves()[valveI].bottomPatchID().active())
        {
            maxLayer = max(maxLayer, valves()[valveI].minBottomLayer());
        }
    
        scalar valveDisplacement = valves_[valveI].curVelocity()*valves_[valveI].cs().axis().z()*engTime().deltaT().value()  ;
        if(valvePistonPosition()[valveI] + engTime().pistonDisplacement().value() >
        valveBottomPosition_[valveI] + valveDisplacement - 5.0*maxLayer - 0.001 )
        {
            deformationValve = true;
        }
    }    

    if(deformationValve)
    {
        return true;
    }
    else  if (virtualPistonPosition() + engTime().pistonDisplacement().value() > deckHeight_ - SMALL)
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
Foam::thoboisSliding::thoboisSliding
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    valves_(*this, engTime().engineDict().lookup("thoboisSliding")),
    movingPointsMaskTopPtr_(NULL),
    movingPointsMaskBottomPtr_(NULL),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
    valveTopTol_(readScalar(engTime().engineDict().lookup("valveTopTol"))),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    valveTopPosition_(nValves(),-GREAT),
    valveBottomPosition_(nValves(),GREAT), 
    valvePistonPosition_(nValves(),GREAT), 
    deckHeight_(GREAT),
    minValveZ_(nValves()),
    poppetValveTol_(readScalar(engTime().engineDict().lookup("poppetValveTol"))),
    bottomValveTol_(readScalar(engTime().engineDict().lookup("bottomValveTol"))),
    msPtr_(motionSolver::New(*this)),
    isReallyClosed_(valves().size(), false),
    correctPointsMotion_(engTime().engineDict().lookup("correctPointsMotion"))

    
{
    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::thoboisSliding::setBoundaryVelocity(volVectorField& U)
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
            U.boundaryField()[valves()[valveI].curtainInPortPatchID().index()] ==
//                valveVel;
                vector::zero;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].curtainInCylinderPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].curtainInCylinderPatchID().index()] ==
//                valveVel;
                vector::zero;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].poppetPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].poppetPatchID().index()] ==
                valveVel;
        }

        if (valves()[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].bottomPatchID().index()] ==
                valveVel;
        }

        
        // If valve is present in geometry, set the motion
        if (valves()[valveI].stemPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].stemPatchID().index()] ==
                valveVel;
        }

     
    }

}

bool Foam::thoboisSliding::inValve(const point& p, const label& i) const
{
    scalar valveX = valves_[i].cs().origin().x();
    scalar valveY = valves_[i].cs().origin().y();
    return (sqrt(sqr(p.x()-valveX)+sqr(p.y()-valveY)) < 0.5*valves_[i].diameter());
}

bool Foam::thoboisSliding::inPiston(const point& p) const
{
    scalar pistonX = piston_.cs().origin().x();
    scalar pistonY = piston_.cs().origin().y();
    return (sqrt(sqr(p.x()-pistonX)+sqr(p.y()-pistonY)) < 0.5*engTime().bore().value());
}


bool Foam::thoboisSliding::isACylinderHeadFace
(
    const labelList& cylHeadFaces, 
    const label face
)
{
    forAll(cylHeadFaces, i)
    {
        if(cylHeadFaces[i] == face)
        {
            return true;
        }
    }
    
    return false;
}



// ************************************************************************* //
