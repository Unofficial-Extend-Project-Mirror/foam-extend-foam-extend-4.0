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
    triSurfaceExtrude2DEdges

Description
    Extrudes edges in x-y plane into a triangulation used for 2D meshing

SourceFiles
    triSurfaceExtrude2DEdges.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceExtrude2DEdges_H
#define triSurfaceExtrude2DEdges_H

#include "triSurf.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceExtrude2DEdges Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceExtrude2DEdges
{
    // Private data
        //- reference to triSurf
        const triSurf& surf_;

    // Private member functions

        //- Disallow default bitwise copy construct
        triSurfaceExtrude2DEdges(const triSurfaceExtrude2DEdges&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceExtrude2DEdges&);

public:

    // Constructors

        //- Construct from octree
        triSurfaceExtrude2DEdges(const triSurf& surface);

    // Destructor

        ~triSurfaceExtrude2DEdges();


    // Member Functions

        //- extrudes edges into a triangulation
        void extrudeSurface(triSurf&) const;

        const triSurf* extrudeSurface() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
