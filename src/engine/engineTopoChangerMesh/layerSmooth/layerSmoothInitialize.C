// The FOAM Project // File: layerSmooth.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   layerSmooth
   \\  /           | Family: engine
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR

-------------------------------------------------------------------------------
*/

#include "layerSmooth.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::layerSmooth::checkAndCalculate()
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
        FatalErrorIn("Foam::layerSmooth::checkAndCalculate()")
            << " : cannot find piston patch"
            << abort(FatalError);
    }

    if (!foundLiner)
    { 
        FatalErrorIn("Foam::layerSmooth::checkAndCalculate()")
            << " : cannot find liner patch"
            << abort(FatalError);
    }

    if (!foundCylinderHead)
    { 
        FatalErrorIn("Foam::layerSmooth::checkAndCalculate()")
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

 /*        
           deckHeight() = max
            (
                boundary()[linerIndex].patch().localPoints()
            ).z();
*/                        
        }
        reduce(deckHeight(), minOp<scalar>());

        Info<< "deckHeight: " << deckHeight() << nl
            << "piston position: " << pistonPosition() << endl;
    }
        
    
} 

void Foam::layerSmooth::setVirtualPositions()
{

    {
        virtualPistonPosition() = -GREAT;

        label pistonFaceIndex = faceZones().findZoneID("pistonLayerFaces");
         
        bool foundPistonFace = (pistonFaceIndex != -1);
        
        if(!foundPistonFace)
        {
            FatalErrorIn("Foam::layerSmooth::setVirtualPistonPosition()")
                << " : cannot find the pistonLayerFaces"
                << exit(FatalError);
    
        }
        
        const labelList& pistonFaces = faceZones()[pistonFaceIndex];
        forAll(pistonFaces, i)
        {
            const face& f = faces()[pistonFaces[i]];
        
            // should loop over facepoints...
            forAll(f, j)
            {
                virtualPistonPosition() = max(virtualPistonPosition(), points()[f[j]].z());
            }
        }
    
        reduce(virtualPistonPosition(), maxOp<scalar>());
    
    }
    

}
