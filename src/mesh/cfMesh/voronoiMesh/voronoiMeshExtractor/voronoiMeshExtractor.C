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
    pointEdgesPtr_(nullptr),
    edgesPtr_(nullptr),
    edgeTetsPtr_(nullptr),
    boundaryEdgePtr_(nullptr)
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
