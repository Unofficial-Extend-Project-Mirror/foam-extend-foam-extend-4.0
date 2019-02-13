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
    triSurfaceCleanupDuplicates

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfaceCleanupDuplicatesDuplicates.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCleanupDuplicates_H
#define triSurfaceCleanupDuplicates_H

#include "triSurf.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class triSurfaceCleanupDuplicates Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceCleanupDuplicates
{
    // Private data
        //- tolerance for points
        const scalar tolerance_;

        //- reference to triSurf
        triSurf& surf_;

        //- reference to octree
        const meshOctree& octree_;

        //- new triangle labels in case some of them is removed
        labelLongList newTriangleLabel_;

        //- check if the operation has already been performed
        bool done_;

    // Private member functions
        //- Check duplicate triangles
        bool checkDuplicateTriangles();

        //- Check duplicate points
        bool mergeDuplicatePoints();

        //- update triangle mapping
        void updateTriangleLabels(const labelLongList&);

        //- Disallow default bitwise copy construct
        triSurfaceCleanupDuplicates(const triSurfaceCleanupDuplicates&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceCleanupDuplicates&);

public:

    // Constructors

        //- Construct from meshOctree and tolerance
        triSurfaceCleanupDuplicates
        (
            const meshOctree& octree,
            const scalar tol=1e-6
        );

    // Destructor

        ~triSurfaceCleanupDuplicates();


    // Member Functions

        //- merge duplicate vertices and triangles
        void mergeIdentities();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
