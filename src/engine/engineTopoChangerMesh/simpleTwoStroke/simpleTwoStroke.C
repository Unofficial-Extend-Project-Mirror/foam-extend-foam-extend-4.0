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
#include "layerAdditionRemoval.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "engineTime.H"
#include "pointField.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleTwoStroke, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, simpleTwoStroke, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



    
bool Foam::simpleTwoStroke::realDeformation() const
{
 

    if (virtualPistonPosition()+ engTime().pistonDisplacement().value() > deckHeight()-engTime().clearance().value()-SMALL)
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
Foam::simpleTwoStroke::simpleTwoStroke
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    scavInCylPatchName_(engTime().engineDict().lookup("scavInCylPatch")),
    scavInPortPatchName_(engTime().engineDict().lookup("scavInPortPatch")),
    movingPointsMaskPtr_(NULL),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
    delta_(readScalar(engTime().engineDict().lookup("delta"))),
    offSet_(readScalar(engTime().engineDict().lookup("offSet"))),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    deckHeight_(GREAT),
    csPtr_
    (
        coordinateSystem::New
        (
            "coordinateSystem",
            engTime().engineDict().subDict("coordinateSystem")
        )
    ),
    foundScavPorts_
    (
        boundaryMesh().findPatchID(scavInCylPatchName_) != -1
        &&
        boundaryMesh().findPatchID(scavInPortPatchName_) != -1
    ),
    scavPortsTol_(readScalar(engTime().engineDict().lookup("scavPortsTol")))
{
    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::simpleTwoStroke::setBoundaryVelocity(volVectorField& U)
{
    vector pistonVel = piston().cs().axis()*engTime().pistonSpeed().value();

    // On the piston movingWallVelocity is used. There is no need to update
    // the piston velocity

//    U.boundaryField()[piston().patchID().index()] = pistonVel;

    Info << "setting the boundary velocity" << endl;
    U.boundaryField()[boundaryMesh().findPatchID(scavInPortPatchName_)] ==
        pistonVel;

}

bool Foam::simpleTwoStroke::portsOpened() const
{

    scalar maxScavPortsZ = max
       (
           boundary()
           [
               boundaryMesh().findPatchID(scavInPortPatchName_)
           ].patch().localPoints()
       ).z();


    scalar minLinerZ = min
       (
           boundary()
           [
               boundaryMesh().findPatchID(scavInCylPatchName_)
           ].patch().localPoints()
       ).z();


    if (maxScavPortsZ - scavPortsTol_ < minLinerZ)
    {
        return false;
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
