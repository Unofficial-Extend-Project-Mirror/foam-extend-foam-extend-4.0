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

Class
    engineTopoChangerMesh

Description
    Base class for engine mesh with topological changes

SourceFiles
    engineTopoChangerMesh.C
    newTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef attachDetachFunctions_H
#define attachDetachFunctions_H

#include "argList.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "faceSet.H"
#include "attachDetach.H"
#include "regionSide.H"
#include "directTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class nearestEqOp
    {

        public:

        void operator()(vector& x, const vector& y) const
        {
            if (magSqr(y) < magSqr(x))
            {
                x = y;
            }
        }
    };


    // Find edge between points v0 and v1.
    label findEdge
    (
       const primitiveMesh& mesh,
       const label v0,
       const label v1
    );


// Checks whether patch present
    void checkPatch
    (
        const polyBoundaryMesh& bMesh,
        const word& name
    );

    void changePatchID
    (
        const polyMesh& mesh,
        const label faceID,
        const label patchID,
        directTopoChange& meshMod
    );



}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
