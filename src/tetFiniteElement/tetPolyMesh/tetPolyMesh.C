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
    The tetrahedral decomposition of the mesh based on the underlying polyMesh.
    Uses minimal storage - the tet decomposition is provided on a cell-by-cell
    basis

\*---------------------------------------------------------------------------*/

#include "tetPolyMesh.H"
#include "demandDrivenData.H"
#include "tetPolyPatchFields.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "tetFemMatrices.H"
#include "mapPolyMesh.H"
#include "MapTetFemFields.H"
#include "tetPolyMeshMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::tetPolyMesh, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void tetPolyMesh::clearOut() const
{
    nPoints_ = -1;
    nEdges_ = -1;
    nTets_ = -1;
    deleteDemandDrivenData(lduPtr_);
    maxNPointsForCell_ = -1;

    clearOutParPointData();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
tetPolyMesh::tetPolyMesh(const polyMesh& pMesh)
:
    GeoMesh<polyMesh>(pMesh),
    tetFemSolution(pMesh),
    boundary_(*this, pMesh.boundaryMesh()),
    faceOffset_(mesh_.nPoints()),
    cellOffset_(faceOffset_ + mesh_.nFaces()),
    nPoints_(-1),
    nEdges_(-1),
    nTets_(-1),
    lduPtr_(NULL),
    maxNPointsForCell_(-1),
    parPointsPtr_(NULL),
    parEdgesPtr_(NULL)
{
    if (debug)
    {
        Info<< "tetPolyMesh::tetPolyMesh(const polyMesh&) : "
            << "Creating tetPolyMesh" << endl;
    }

    addParallelPointPatch();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetPolyMesh::~tetPolyMesh()
{
    if (debug)
    {
        Info<< "tetPolyMesh::~tetPolyMesh() : "
            << "Deleting tetPolyMesh" << endl;
    }

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return number of edges in decomposition for a face
label tetPolyMesh::nEdgesForFace(const label faceID) const
{
    // It is all the edges around the circumference + the internal diags

    // Note: compiler bug.  It should be possible to do this without wrapping
    const faceList& f = mesh_.faces();
    return 2*f[faceID].size();
}


// Return number of edges in decomposition connected to a given point
label tetPolyMesh::nEdgesForPoint(const label pointID) const
{
    const label startFaceOwn = lduAddr().ownerStartAddr()[pointID];
    const label endFaceOwn = lduAddr().ownerStartAddr()[pointID + 1];

    const label startFaceNbr = lduAddr().losortStartAddr()[pointID];
    const label endFaceNbr = lduAddr().losortStartAddr()[pointID + 1];

    // pointID is the owner of the first lot and the neighbour of the second lot
    return (endFaceOwn - startFaceOwn) + (endFaceNbr - startFaceNbr);
}


// Return number of edges in decomposition connected to a given point
labelList tetPolyMesh::edgesForPoint(const label pointID) const
{
    const label startFaceOwn = lduAddr().ownerStartAddr()[pointID];
    const label endFaceOwn = lduAddr().ownerStartAddr()[pointID + 1];

    const label startFaceNbr = lduAddr().losortStartAddr()[pointID];
    const label endFaceNbr = lduAddr().losortStartAddr()[pointID + 1];

    const unallocLabelList& losort = lduAddr().losortAddr();

    // pointID is the owner of the first lot and the neighbour of the second lot
    labelList edgeIndices(nEdgesForPoint(pointID), -1);

    label i = 0;

    // owner side
    for (label edgeLabel = startFaceOwn; edgeLabel < endFaceOwn; edgeLabel++)
    {
        edgeIndices[i] = edgeLabel;
        i++;
    }

    // neighbour side
    for (label edgeLabel = startFaceNbr; edgeLabel < endFaceNbr; edgeLabel++)
    {
        edgeIndices[i] = losort[edgeLabel];
        i++;
    }

    return edgeIndices;
}


void tetPolyMesh::updateMesh
(
    const tetPolyMeshMapper& mapper
)
{
    // It is assumed that polyMesh::morph() has already been called.
    if (debug)
    {
        Info<< "void tetPolyMesh::updateMesh: "
            << "Mesh update on topological change" << endl;
    }

    // Clear out all data
    clearOut();

    // Reset face and cell offset (topology has changed)
    faceOffset_ = mesh_.nPoints();
    cellOffset_ = faceOffset_ + mesh_.nFaces();

    boundary_.updateMesh();

    // Map all the tetPointFields in the database

    MapGeometricFields
    <
        scalar,
        tetPolyPatchField,
        tetPolyMeshMapper,
        tetPointMesh
    >(mapper);

    MapGeometricFields
    <
        vector,
        tetPolyPatchField,
        tetPolyMeshMapper,
        tetPointMesh
    >(mapper);

    MapGeometricFields
    <
        tensor,
        tetPolyPatchField,
        tetPolyMeshMapper,
        tetPointMesh
    >(mapper);

    // Map all the element fields in the database

    MapGeometricFields
    <
        scalar,
        elementPatchField,
        tetPolyMeshMapper,
        elementMesh
    >(mapper);

    MapGeometricFields
    <
        vector,
        elementPatchField,
        tetPolyMeshMapper,
        elementMesh
    >(mapper);

    MapGeometricFields
    <
        tensor,
        elementPatchField,
        tetPolyMeshMapper,
        elementMesh
    >(mapper);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool tetPolyMesh::operator!=(const tetPolyMesh& bm) const
{
    return &bm != this;
}


bool tetPolyMesh::operator==(const tetPolyMesh& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
