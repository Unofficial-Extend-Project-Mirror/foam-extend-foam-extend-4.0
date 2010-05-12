/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Description
    The tetrahedral decomposition of the mesh based on the underlying polyMesh.
    Uses minimal storage - the tet decomposition is provided on a cell-by-cell
    basis

\*---------------------------------------------------------------------------*/

#include "tetPolyMeshFaceDecomp.H"
#include "demandDrivenData.H"
#include "tetPolyPatchFields.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "tetFemMatrices.H"
#include "mapPolyMesh.H"
#include "MapTetFemFields.H"
#include "tetPolyMeshMapperFaceDecomp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::tetPolyMeshFaceDecomp, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void tetPolyMeshFaceDecomp::clearOut() const
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
tetPolyMeshFaceDecomp::tetPolyMeshFaceDecomp(const polyMesh& pMesh)
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

tetPolyMeshFaceDecomp::~tetPolyMeshFaceDecomp()
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
label tetPolyMeshFaceDecomp::nEdgesForFace(const label faceID) const
{
    // It is all the edges around the circumference + the internal diags

    // Note: compiler bug.  It should be possible to do this without wrapping
    const faceList& f = mesh_.faces();
    return 2*f[faceID].size();
}


// Return number of edges in decomposition connected to a given point
label tetPolyMeshFaceDecomp::nEdgesForPoint(const label pointID) const
{
    const label startFaceOwn = lduAddr().ownerStartAddr()[pointID];
    const label endFaceOwn = lduAddr().ownerStartAddr()[pointID + 1];

    const label startFaceNbr = lduAddr().losortStartAddr()[pointID];
    const label endFaceNbr = lduAddr().losortStartAddr()[pointID + 1];

    // pointID is the owner of the first lot and the neighbour of the second lot
    return (endFaceOwn - startFaceOwn) + (endFaceNbr - startFaceNbr);
}


// Return number of edges in decomposition connected to a given point
labelList tetPolyMeshFaceDecomp::edgesForPoint(const label pointID) const
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


void tetPolyMeshFaceDecomp::updateMesh
(
    const tetPolyMeshMapperFaceDecomp& mapper
)
{
    // It is assumed that polyMesh::morph() has already been called.
    if (debug)
    {
        Info<< "void tetPolyMeshFaceDecomp::updateMesh: "
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
        tetPolyMeshMapperFaceDecomp,
        tetPointMesh
    >(mapper);

    MapGeometricFields
    <
        vector,
        tetPolyPatchField,
        tetPolyMeshMapperFaceDecomp,
        tetPointMesh
    >(mapper);

    MapGeometricFields
    <
        tensor,
        tetPolyPatchField,
        tetPolyMeshMapperFaceDecomp,
        tetPointMesh
    >(mapper);

    // Map all the element fields in the database

    MapGeometricFields
    <
        scalar,
        elementPatchField,
        tetPolyMeshMapperFaceDecomp,
        elementMesh
    >(mapper);

    MapGeometricFields
    <
        vector,
        elementPatchField,
        tetPolyMeshMapperFaceDecomp,
        elementMesh
    >(mapper);

    MapGeometricFields
    <
        tensor,
        elementPatchField,
        tetPolyMeshMapperFaceDecomp,
        elementMesh
    >(mapper);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool tetPolyMeshFaceDecomp::operator!=(const tetPolyMeshFaceDecomp& bm) const
{
    return &bm != this;
}


bool tetPolyMeshFaceDecomp::operator==(const tetPolyMeshFaceDecomp& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
