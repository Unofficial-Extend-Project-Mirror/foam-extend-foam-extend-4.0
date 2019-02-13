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

Class
    findCellsIntersectingSurface

Description
    Finds which cells in the mesh intersect the surface mesh

SourceFiles
    findCellsIntersectingSurface.C

\*---------------------------------------------------------------------------*/

#ifndef findCellsIntersectingSurface_H
#define findCellsIntersectingSurface_H

#include "boolList.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMeshGen;
class triSurf;
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class findCellsIntersectingSurface Declaration
\*---------------------------------------------------------------------------*/

class findCellsIntersectingSurface
{
    // Private data
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- Pointer to the octree
        meshOctree* octreePtr_;

        //- check whether the octree was generated or not
        const bool octreeGenerated_;

        //- stores information about intersected cells
        boolList intersectedCells_;

        //- stores information which surface facets intersect each cell
        VRWGraph facetsIntersectingCell_;

    // Private member functions
        //- generate the octree
        void generateOctree(const triSurf&);

        //- check for the intersected cells
        void findIntersectedCells();

public:

    // Constructors

        //- Construct from mesh and octree
        findCellsIntersectingSurface
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );

        //- Construct from mesh and surface
        findCellsIntersectingSurface
        (
            polyMeshGen& mesh,
            const triSurf& surface
        );

    // Destructor
        ~findCellsIntersectingSurface();

    // Public member functions
        //- return the list of intersected cells;
        const boolList& intersectedCells() const;

        //- return the graph of facets intersecting each cell
        const VRWGraph& facetsIntersectingCells() const;

        //- create a cell subset containing intersected cells
        void addIntersectedCellsToSubset
        (
            const word subsetName = "intersectedCells"
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
