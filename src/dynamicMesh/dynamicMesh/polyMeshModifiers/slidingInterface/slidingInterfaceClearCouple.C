/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Class
    slidingInterface

Description
    Sliding interface clear-out

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.  Copyright Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#include "slidingInterface.H"
#include "polyTopoChange.H"
#include "polyMesh.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::clearCouple
(
    polyTopoChange& ref
) const
{
    if (debug)
    {
        Pout<< "void slidingInterface::clearCouple("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Clearing old couple points and faces." << endl;
    }

    // Remove all points from the point zone

    const polyMesh& mesh = topoChanger().mesh();

    const labelList& cutPointZoneLabels =
        mesh.pointZones()[cutPointZoneID_.index()];

    forAll (cutPointZoneLabels, pointI)
    {
//         Pout << "Removing point in decouple " << cutPointZoneLabels[pointI] << endl;
        ref.setAction(polyRemovePoint(cutPointZoneLabels[pointI]));
    }

    // Remove all faces from the face zone
    const labelList& cutFaceZoneLabels =
        mesh.faceZones()[cutFaceZoneID_.index()];

    forAll (cutFaceZoneLabels, faceI)
    {
        ref.setAction(polyRemoveFace(cutFaceZoneLabels[faceI]));
    }

    if (debug)
    {
        Pout<< "void slidingInterface::clearCouple("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Finished clearing old couple points and faces." << endl;
    }
}


// ************************************************************************* //
