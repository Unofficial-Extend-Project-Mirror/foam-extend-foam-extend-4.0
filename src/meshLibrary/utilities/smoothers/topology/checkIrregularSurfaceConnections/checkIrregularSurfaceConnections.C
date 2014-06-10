/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "checkIrregularSurfaceConnections.H"

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkIrregularSurfaceConnections::checkIrregularSurfaceConnections
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    meshSurfacePtr_(NULL)
{
}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkIrregularSurfaceConnections::~checkIrregularSurfaceConnections()
{
    clearMeshEngine();
    
    mesh_.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkIrregularSurfaceConnections::checkIrregularVertices
(
    labelHashSet& badVertices
)
{
    checkAndFixCellGroupsAtBndVertices(badVertices, false);
    
    checkEdgeFaceConnections(badVertices, false);
        
    checkFaceGroupsAtBndVertices(badVertices, false);
}

bool checkIrregularSurfaceConnections::checkAndFixIrregularConnections()
{
    Info << "Checking for irregular surface connections" << endl;
    
    bool finished;
    
    labelHashSet badVertices;
    
    do
    {
        finished = true;
        
        while( checkAndFixCellGroupsAtBndVertices(badVertices, true) )
            finished = false;
        
        while( checkEdgeFaceConnections(badVertices, true) )
            finished = false;
        
        if( checkFaceGroupsAtBndVertices(badVertices, true) )
            finished = false;
    } while( !finished );
    
    polyMeshGenModifier(mesh_).removeUnusedVertices();
    
    Info << "Finished checking for irregular surface connections" << endl;
    
    if( returnReduce(badVertices.size(), sumOp<label>()) != 0 )
        return true;
    
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
