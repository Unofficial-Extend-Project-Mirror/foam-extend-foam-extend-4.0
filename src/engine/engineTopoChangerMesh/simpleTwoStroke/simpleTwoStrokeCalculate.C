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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleTwoStroke::correctVerticalMotion()
{
}


void Foam::simpleTwoStroke::calcMovingMasks() const
{
    if (debug)
    {
        Info<< "void movingSquaresTM::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void movingSquaresTM::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    const labelList& cellAddr =
        cellZones()[cellZones().findZoneID("movingCells")];

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    
    if(foundScavPorts())
    {

        const word innerScavZoneName
        (
            scavInCylPatchName_  + "Zone"
        );

        const labelList& innerScavAddr =
            faceZones()[faceZones().findZoneID(innerScavZoneName)];

        forAll (innerScavAddr, faceI)
        {
            const face& curFace = f[innerScavAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }

        const word outerScavZoneName
        (
            scavInPortPatchName_ + "Zone"
        );

        const labelList& outerScavAddr =
            faceZones()[faceZones().findZoneID(outerScavZoneName)];

        forAll (outerScavAddr, faceI)
        {
            const face& curFace = f[outerScavAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 0;
            }
        }
    
    }
           
}

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::simpleTwoStroke::movingPointsMask() const
{
    if(movingPointsMaskPtr_)
    {
        movingPointsMaskPtr_ = NULL;
    }
    
    if (!movingPointsMaskPtr_)
    {
        calcMovingMasks();
    }

    return *movingPointsMaskPtr_;
}

void Foam::simpleTwoStroke::applyMovingMask()
{
}

