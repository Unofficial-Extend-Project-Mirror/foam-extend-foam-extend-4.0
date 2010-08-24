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

#include "pistonSliding.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "componentMixedTetPolyPatchVectorField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pistonSliding::calcMovingMaskTop(const label i) const
{
    if (debug)
    {
        Info<< "void movingSquaresTM::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskTopPtr_)
    {
        FatalErrorIn("void movingSquaresTM::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskTopPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMaskTop = *movingPointsMaskTopPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    const labelList& cellTopVAddr =
        cellZones()[cellZones().findZoneID("movingCellsTopZoneV"+ Foam::name(i+1))];
    
    bool foundOne(false);

    forAll (cellTopVAddr, cellI)
    {
        const cell& curCell = c[cellTopVAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                

                {
                    foundOne = true;  
                    movingPointsMaskTop[curFace[pointI]] = 1;
                }
                
                
            }
            
        }
        
    }
    
    const labelList& cellTopAddr =
        cellZones()[cellZones().findZoneID("movingCellsZoneV"+ Foam::name(i+1))];
    
    forAll (cellTopAddr, cellI)
    {
        const cell& curCell = c[cellTopAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                

                {
                    foundOne = true;  
                    movingPointsMaskTop[curFace[pointI]] = 1;
                }
                
                
            }
            
        }
        
    }
    
    
    
/*    
    if(valves_[i].poppetPatchID().active())
    {

        const word curtainCylZoneName
        (
            "curtainCylZoneV" + Foam::name(i + 1)
        );

        const labelList& curtainCylAddr =
            faceZones()[faceZones().findZoneID(curtainCylZoneName)];

        forAll (curtainCylAddr, faceI)
        {
            const face& curFace = f[curtainCylAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMaskTop[curFace[pointI]] = 0;
            }
        }
    
    }
*/    
}

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::pistonSliding::movingPointsMaskTop(const label i) const
{
    if(movingPointsMaskTopPtr_)
    {
        movingPointsMaskTopPtr_ = NULL;
    }
    
    if (!movingPointsMaskTopPtr_)
    {
        calcMovingMaskTop(i);
    }

    return *movingPointsMaskTopPtr_;
}


void Foam::pistonSliding::calcMovingMaskBottom(const label i) const
{
    if (debug)
    {
        Info<< "void movingSquaresTM::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskBottomPtr_)
    {
        FatalErrorIn("void movingSquaresTM::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskBottomPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMaskBottom = *movingPointsMaskBottomPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    const labelList& cellAddr =
        cellZones()[cellZones().findZoneID("movingCellsBottomZoneV"+ Foam::name(i+1))];

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                {
                    movingPointsMaskBottom[curFace[pointI]] = 1;
                }
            }
        }
    }

/*    
    if(valves_[i].bottomPatchID().active())
    {

        const word curtainCylZoneName
        (
            "curtainCylZoneV" + Foam::name(i + 1)
        );

        const labelList& curtainCylAddr =
            faceZones()[faceZones().findZoneID(curtainCylZoneName)];

        forAll (curtainCylAddr, faceI)
        {
            const face& curFace = f[curtainCylAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMaskBottom[curFace[pointI]] = 0;
            }
        }
    
    }
*/    
}

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::pistonSliding::movingPointsMaskBottom(const label i) const
{
    if(movingPointsMaskBottomPtr_)
    {
        movingPointsMaskBottomPtr_ = NULL;
    }
    
    if (!movingPointsMaskBottomPtr_)
    {
        calcMovingMaskBottom(i);
    }

    return *movingPointsMaskBottomPtr_;
}


void Foam::pistonSliding::calcMovingMaskPiston() const
{
    if (debug)
    {
        Info<< "void movingSquaresTM::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskPistonPtr_)
    {
        FatalErrorIn("void movingSquaresTM::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPistonPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMaskPiston = *movingPointsMaskPistonPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    const labelList& cellAddr =
        cellZones()[cellZones().findZoneID("movingCellsPiston")];
    
    bool foundOne(false);

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                

                {
                    foundOne = true;  
                    movingPointsMaskPiston[curFace[pointI]] = 1;
                }
                
                
            }
            
        }
        
    }
    
        
}

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::pistonSliding::movingPointsMaskPiston() const
{
    if(movingPointsMaskPistonPtr_)
    {
        movingPointsMaskPistonPtr_ = NULL;
    }
    
    if (!movingPointsMaskPistonPtr_)
    {
        calcMovingMaskPiston();
    }

    return *movingPointsMaskPistonPtr_;
}



void Foam::pistonSliding::calcMovingMaskPistonValves(const label i) const
{
    if (debug)
    {
        Info<< "void movingSquaresTM::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskPistonValvesPtr_)
    {
        FatalErrorIn("void movingSquaresTM::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPistonValvesPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMaskPistonValves = *movingPointsMaskPistonValvesPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    const labelList& cellAddr =
        cellZones()[cellZones().findZoneID("movingCellsPistonV" + Foam::name(i+1))];
    
    bool foundOne(false);

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                

                {
                    foundOne = true;  
                    movingPointsMaskPistonValves[curFace[pointI]] = 1;
                }
                
                
            }
            
        }
        
    }

        
}

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::pistonSliding::movingPointsMaskPistonValves(const label i) const
{
    if(movingPointsMaskPistonValvesPtr_)
    {
        movingPointsMaskPistonValvesPtr_ = NULL;
    }
    
    if (!movingPointsMaskPistonValvesPtr_)
    {
        calcMovingMaskPistonValves(i);
    }

    return *movingPointsMaskPistonValvesPtr_;
}
