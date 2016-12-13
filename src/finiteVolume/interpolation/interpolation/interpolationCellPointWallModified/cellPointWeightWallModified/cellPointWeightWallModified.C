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

\*---------------------------------------------------------------------------*/

#include "cellPointWeightWallModified.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointWeightWallModified::cellPointWeightWallModified
(
    const polyMesh& mesh,
    const vector& position,
    const label cellIndex,
    const label faceIndex
)
:
    cellPointWeight(mesh, position, cellIndex, faceIndex)
{
    if (faceIndex < 0)
    {
        findTetrahedron(mesh, position, cellIndex);
    }
    else
    {
        const polyBoundaryMesh& bm = mesh.boundaryMesh();
        label patchI = bm.whichPatch(faceIndex);
        if (patchI != -1)
        {
            if (bm[patchI].isWall())
            {
                // Apply cell centre value wall faces
                weights_[0] = 0.0;
                weights_[1] = 0.0;
                weights_[2] = 0.0;
                weights_[3] = 1.0;
            }
        }
        else
        {
            // Interpolate
            findTriangle(mesh, position, faceIndex);
        }
    }
}


// ************************************************************************* //
