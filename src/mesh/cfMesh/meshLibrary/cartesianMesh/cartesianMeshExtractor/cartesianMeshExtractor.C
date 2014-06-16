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

#include "cartesianMeshExtractor.H"
#include "meshOctree.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshExtractor::clearOut()
{
    deleteDemandDrivenData(leafCellLabelPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
cartesianMeshExtractor::cartesianMeshExtractor
(
    meshOctree& octree,
    const IOdictionary& meshDict,
    polyMeshGen& mesh
)
:
    octreeCheck_(octree, meshDict, false),
    mesh_(mesh),
    decomposeSplitHexes_(false),
    leafCellLabelPtr_(new labelList(octree.numberOfLeaves(), -1))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cartesianMeshExtractor::~cartesianMeshExtractor()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshExtractor::decomposeSplitHexes()
{
    decomposeSplitHexes_ = true;
}

void cartesianMeshExtractor::createMesh()
{
    Info << "Extracting polyMesh" << endl;
    
    //- create points and pointLeaves addressing
    createPointsAndAddressing();

    //- create the mesh
    createPolyMesh();
    
    //- decompose split-hex cells into tetrahedra and pyramids
    decomposeSplitHexesIntoTetsAndPyramids();
    
    //- remove unused vertices
    polyMeshGenModifier(mesh_).removeUnusedVertices();
    
    Info << "Mesh has :" << nl
        << mesh_.points().size() << " vertices " << nl
        << mesh_.faces().size() << " faces" << nl
        << mesh_.cells().size() << " cells" << endl;
    
    if( Pstream::parRun() )
    {
        label nCells = mesh_.cells().size();
        reduce(nCells, sumOp<label>());
        Info << "Total number of cells " << nCells << endl;
    }
    if( mesh_.cells().size() == 0 )
    {
        FatalErrorIn
        (
            "void cartesianMeshExtractor::createMesh()"
        ) << "There are no cells in the mesh!"
        << nl << "The reasons for this can be fwofold:"
        << nl << "1. Inadequate mesh resolution."
        << nl << "2. You maxCellSize is a multiplier of the domain length."
        << " This can be reolved by reducing the maxCellSize by a fraction."
        << "i.e. 2.49999 instead of 2.5." << exit(FatalError);
    }
    
    Info << "Finished extracting polyMesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
