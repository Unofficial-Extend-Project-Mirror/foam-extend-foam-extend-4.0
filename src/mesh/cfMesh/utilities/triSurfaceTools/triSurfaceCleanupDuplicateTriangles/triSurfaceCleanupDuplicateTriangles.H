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
    Removes duplicate triangles from the surface.

SourceFiles
    triSurfaceCleanupDuplicateTriangles.C
    triSurfaceCleanupDuplicateTrianglesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCleanupDuplicateTriangles_H
#define triSurfaceCleanupDuplicateTriangles_H

#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;

/*---------------------------------------------------------------------------*\
            Class triSurfaceCleanupDuplicateTriangles Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceCleanupDuplicateTriangles
{
    // Private data
        //- reference to triSurf
        triSurf& surf_;

        //- new triangle labels in case some of them is removed
        labelLongList newTriangleLabel_;

    // Private member functions
        //- Check duplicate triangles
        void checkDuplicateTriangles();

        //- Disallow default bitwise copy construct
        triSurfaceCleanupDuplicateTriangles
        (
            const triSurfaceCleanupDuplicateTriangles&
        );

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceCleanupDuplicateTriangles&);

public:

    // Constructors

        //- Construct from triSurf
        triSurfaceCleanupDuplicateTriangles(triSurf&);

    // Destructor

        ~triSurfaceCleanupDuplicateTriangles();


    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
