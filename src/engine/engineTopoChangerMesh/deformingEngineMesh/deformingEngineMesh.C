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

#include "deformingEngineMesh.H"
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
#include "tetMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deformingEngineMesh, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, deformingEngineMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::deformingEngineMesh::realDeformation() const
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
Foam::deformingEngineMesh::deformingEngineMesh
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    valves_(*this, engTime().engineDict().lookup("deformingEngineMesh")),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    deckHeight_(GREAT),
    cylinderHeadName_(engTime().engineDict().lookup("cylinderHeadName")),
    linerName_(engTime().engineDict().lookup("linerName")),
    headCellSetName_(engTime().engineDict().lookup("headCellSetName"))

{
    // Add zones and modifiers if not already there.
    addMeshZones();
    msPtr_ = motionSolver::New(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::deformingEngineMesh::setBoundaryVelocity(volVectorField& U)
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

    }

}

// ************************************************************************* //
