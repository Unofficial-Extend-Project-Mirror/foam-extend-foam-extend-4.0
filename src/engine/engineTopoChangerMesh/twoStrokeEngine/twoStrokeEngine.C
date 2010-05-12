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
    defineTypeNameAndDebug(twoStrokeEngine, 0);
    addToRunTimeSelectionTable
    (
        engineTopoChangerMesh,
        twoStrokeEngine,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::twoStrokeEngine::realDeformation() const
{
    if
    (
        virtualPistonPosition() + engTime().pistonDisplacement().value()
      > deckHeight() - engTime().clearance().value() - SMALL
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
Foam::twoStrokeEngine::twoStrokeEngine
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    scavInPortPatches_(engTime().engineDict().lookup("scavInPortPatches")),
    scavInCylPatches_(engTime().engineDict().lookup("scavInCylPatches")),
    headPointsSetName_(engTime().engineDict().lookup("headPointsSetName")),
    headCellsSetName_(engTime().engineDict().lookup("headCellsSetName")),
    movingCellSetName_(engTime().engineDict().lookup("movingCellSetName")),
    movingPointsMaskPtr_(NULL),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
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
    )

{
    if(scavInPortPatches_.size() != scavInCylPatches_.size())
    {
            FatalErrorIn
            (
                "Foam::twoStrokeEngine::twoStrokeEngine(const IOobject& io)"
            )   << "The size of the scavenging-cylinder patches is not"
                << "the same of the scavenging-port patches"
                << abort(FatalError);
    }

    forAll(scavInPortPatches_, patchi)
    {
        if(boundaryMesh().findPatchID(scavInPortPatches_[patchi]) == -1)
        {
            FatalErrorIn
            (
                "Foam::twoStrokeEngine::twoStrokeEngine(const IOobject& io)"
            )   << "patch called" << scavInPortPatches_[patchi]
                << "does not exist"
                << abort(FatalError);
        }
    }

    forAll(scavInCylPatches_, patchi)
    {
        if(boundaryMesh().findPatchID(scavInCylPatches_[patchi]) == -1)
        {
            FatalErrorIn
            (
                "Foam::twoStrokeEngine::twoStrokeEngine(const IOobject& io)"
            )   << "patch called" << scavInCylPatches_[patchi]
                << "does not exist"
                << abort(FatalError);
        }
    }

    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::twoStrokeEngine::setBoundaryVelocity(volVectorField& U)
{
    vector pistonVel = piston().cs().axis()*engTime().pistonSpeed().value();

    //  On the piston movingWallVelocity is used.
    // There is no need to update the piston velocity
//    U.boundaryField()[piston().patchID().index()] = pistonVel;

    forAll(scavInPortPatches_, patchi)
    {
        U.boundaryField()
            [boundaryMesh().findPatchID(scavInPortPatches_[patchi])] ==
            pistonVel;
    }
}


// ************************************************************************* //
