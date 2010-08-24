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

#include "pistonDeform.H"
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
    defineTypeNameAndDebug(pistonDeform, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, pistonDeform, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::pistonDeform::checkAndCalculate()
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
        FatalErrorIn("Foam::pistonDeform::checkAndCalculate()")
            << " : cannot find piston patch"
            << abort(FatalError);
    }

    if (!foundLiner)
    { 
        FatalErrorIn("Foam::pistonDeform::checkAndCalculate()")
            << " : cannot find liner patch"
            << abort(FatalError);
    }

    if (!foundCylinderHead)
    { 
        FatalErrorIn("Foam::pistonDeform::checkAndCalculate()")
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pistonDeform::pistonDeform
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    pistonPosition_(-GREAT),
    deckHeight_(GREAT)
{
    checkAndCalculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pistonDeform::~pistonDeform()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pistonDeform::update()
{
    scalar deltaZ = engTime().pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << endl;
    
    Info << "pistonPosition = " << pistonPosition_ << endl;
    Info << "deckHeight = " << deckHeight_ << endl;
    

    // Position of the top of the static mesh layers above the piston
    scalar pistonPlusLayers = pistonPosition_; //+ pistonLayers_.value();

    pointField newPoints = points();

    forAll (newPoints, pointi)
    {
        point& p = newPoints[pointi];

        if (p.z() < pistonPlusLayers)           // In piston bowl
        {
            p.z() += deltaZ;
        }
        else if (p.z() < deckHeight_)   // In liner region
        {
            p.z() += 
                deltaZ
               *(deckHeight_ - p.z())
               /(deckHeight_ - pistonPlusLayers);
        }
    }

    {
        movePoints(newPoints);
    }

    pistonPosition_ += deltaZ;
    scalar pistonSpeed = deltaZ/engTime().deltaT().value();

    Info<< "clearance: " << deckHeight_ - pistonPosition_ << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;


    return false;

}

void Foam::pistonDeform::setBoundaryVelocity(volVectorField& U)
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
