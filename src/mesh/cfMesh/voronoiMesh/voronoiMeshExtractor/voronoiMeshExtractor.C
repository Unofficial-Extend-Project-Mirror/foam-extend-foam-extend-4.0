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

#include "voronoiMeshExtractor.H"
#include "demandDrivenData.H"

// #define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label voronoiMeshExtractor::sameOrientation_[6] = {3, 1, 2, 2, 3, 0};

label voronoiMeshExtractor::oppositeOrientation_[6] = {2, 3, 1, 0, 0, 1};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::clearOut()
{
    deleteDemandDrivenData(pointEdgesPtr_);
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(edgeTetsPtr_);
    deleteDemandDrivenData(boundaryEdgePtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
voronoiMeshExtractor::voronoiMeshExtractor
(
    const meshOctree& octree,
    const IOdictionary& meshDict,
    polyMeshGen& mesh
)
:
    tetCreator_(octree, meshDict),
    mesh_(mesh),
    pointEdgesPtr_(NULL),
    edgesPtr_(NULL),
    edgeTetsPtr_(NULL),
    boundaryEdgePtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voronoiMeshExtractor::~voronoiMeshExtractor()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::createMesh()
{
    Info << "Extracting voronoi mesh" << endl;

    //- copy tet points into the mesh
    createPoints();

    //- create the mesh
    createPolyMesh();

    polyMeshGenModifier(mesh_).reorderBoundaryFaces();
    polyMeshGenModifier(mesh_).removeUnusedVertices();

    Info << "Mesh has :" << nl
        << mesh_.points().size() << " vertices " << nl
        << mesh_.faces().size() << " faces" << nl
        << mesh_.cells().size() << " cells" << endl;

    Info << "Finished extracting voronoi mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
