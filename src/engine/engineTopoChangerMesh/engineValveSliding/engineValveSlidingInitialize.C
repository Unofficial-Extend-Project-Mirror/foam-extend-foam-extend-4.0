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

#include "engineValveSliding.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::engineValveSliding::checkAndCalculate()
{
    
    label pistonIndex = -1;
    bool foundPiston = false;

    label linerIndex = -1;
    bool foundLiner = false;

    label cylinderHeadIndex = -1;
    bool foundCylinderHead = false;
    
    
    forAll(boundary(), i)
    {
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
        FatalErrorIn("Foam::verticalValvesGambit::checkAndCalculate()")
            << " : cannot find piston patch"
            << abort(FatalError);
    }

    if (!foundLiner)
    { 
        FatalErrorIn("Foam::verticalValvesGambit::checkAndCalculate()")
            << " : cannot find liner patch"
            << abort(FatalError);
    }

    if (!foundCylinderHead)
    { 
        FatalErrorIn("Foam::verticalValvesGambit::checkAndCalculate()")
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

void Foam::engineValveSliding::setVirtualPositions()
{}
