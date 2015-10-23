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

Description

\*---------------------------------------------------------------------------*/

#include "contactPatchPair.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void contactPatchPair::updateContact
(
    const volVectorField& disp
)
{
    vectorField slaveVertexResult =
        patchToPatchInterpolate_.pointInterpolate<vector>
        (
            masterInterpolate_.faceToPointInterpolate
            (
                disp.boundaryField()[masterPatchIndex_]
            )
        );

    const vectorField& projectionDirection =
        mesh_.boundaryMesh()[slavePatchIndex_].pointNormals();

    // use this for touch fraction
    scalarField vertexSlaveGap =
    (
        (
            slaveVertexResult
          - slaveInterpolate_.faceToPointInterpolate
            (
                disp.boundaryField()[slavePatchIndex_]
            )
        )
        & projectionDirection
    ) + patchToPatchInterpolate_.pointDistanceToIntersection() - tol_;

    // use this for slaveDisplacement
    vectorField vertexSlaveDisplacement =
        slaveVertexResult
      + patchToPatchInterpolate_.pointDistanceToIntersection()
        *projectionDirection;

    // calculate area in contact

    const faceList& slavePatchLocalFaces =
        mesh_.boundaryMesh()[slavePatchIndex_].localFaces();

    const pointField& slavePatchLocalPoints =
        mesh_.boundaryMesh()[slavePatchIndex_].localPoints();

    touchFraction_ = 0;

    forAll (slavePatchLocalFaces, faceI)
    {
        touchFraction_[faceI] =
            slavePatchLocalFaces[faceI].areaInContact
            (
                slavePatchLocalPoints,
                vertexSlaveGap
            );
    }

    // calculate contact direction
    slaveDisplacement_ =
        slaveInterpolate_.pointToFaceInterpolate(vertexSlaveDisplacement);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
