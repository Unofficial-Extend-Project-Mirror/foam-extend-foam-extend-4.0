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
    meshSurfaceEdgeExtractorNonTopo

Description
    Stores boundary faces into patches and captures edges and corners

SourceFiles
    meshSurfaceEdgeExtractorNonTopo.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractorNonTopo_H
#define meshSurfaceEdgeExtractorNonTopo_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractorNonTopo Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractorNonTopo
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;

        //- octree
        const meshOctree& meshOctree_;

    // Private member functions

        //- distribute boundary faces into patches
        void distributeBoundaryFaces();

        //- re-map points after edges have been extracted
        void remapBoundaryPoints();

        //- Disallow default bitwise copy construct
        meshSurfaceEdgeExtractorNonTopo(const meshSurfaceEdgeExtractorNonTopo&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractorNonTopo&);

public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractorNonTopo
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );

    // Destructor

        ~meshSurfaceEdgeExtractorNonTopo();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
