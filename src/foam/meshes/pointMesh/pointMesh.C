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

#include "pointMesh.H"
#include "globalMeshData.H"
#include "globalPointPatch.H"
#include "pointMeshMapper.H"
#include "pointFields.H"
#include "MapGeometricFields.H"
#include "MapPointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::pointMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMesh::mapFields(const mapPolyMesh& mpm) const
{
    // Create a mapper
    const pointMeshMapper m(*this, mpm);

    MapGeometricFields<scalar, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields<vector, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields
    <
        sphericalTensor,
        pointPatchField,
        pointMeshMapper,
        pointMesh
    >(m);

    MapGeometricFields<symmTensor, pointPatchField, pointMeshMapper, pointMesh>
    (m);

    MapGeometricFields
    <
        symmTensor4thOrder,
        pointPatchField,
        pointMeshMapper,
        pointMesh
    >(m);

    MapGeometricFields<diagTensor, pointPatchField, pointMeshMapper, pointMesh>
    (m);

    MapGeometricFields<tensor, pointPatchField, pointMeshMapper, pointMesh>(m);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh
(
    const polyMesh& pMesh,
    bool alwaysConstructGlobalPatch
)
:
    MeshObject<polyMesh, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    boundary_(*this, pMesh.boundaryMesh())
{
    // Add the globalPointPatch if there are global points
    if
    (
        alwaysConstructGlobalPatch
     || GeoMesh<polyMesh>::mesh_.globalData().nGlobalPoints()
    )
    {
        boundary_.setSize(boundary_.size() + 1);

        boundary_.set
        (
            boundary_.size() - 1,
            new globalPointPatch
            (
                boundary_,
                boundary_.size() - 1
            )
        );
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


bool Foam::pointMesh::movePoints() const
{
    // Casting const-ness to answer the interface of meshObject
    // HJ, 30/Aug/2010
    const_cast<pointBoundaryMesh&>(boundary_).movePoints();

    return true;
}


bool Foam::pointMesh::updateMesh(const mapPolyMesh& mpm) const
{
    // Casting const-ness to answer the interface of meshObject
    // HJ, 30/Aug/2010
    const_cast<pointBoundaryMesh&>(boundary_).updateMesh();

    // Map all registered point fields
    mapFields(mpm);

    return true;
}


// ************************************************************************* //
