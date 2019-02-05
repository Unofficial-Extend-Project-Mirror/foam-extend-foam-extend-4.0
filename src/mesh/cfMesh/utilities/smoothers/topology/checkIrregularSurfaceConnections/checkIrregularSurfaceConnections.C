/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
    meshSurfacePtr_(nullptr)
{}

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
