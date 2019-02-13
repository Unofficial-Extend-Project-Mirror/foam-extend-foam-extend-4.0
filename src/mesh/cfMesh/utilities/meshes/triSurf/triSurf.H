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
    triSurf

Description
    A class for triangulated surface used in the meshing process.
    It is derived from points and facets with some additional subsets.
    Subsets are vert useful for local mesh refinement.

SourceFiles
    triSurf.C

\*---------------------------------------------------------------------------*/

#ifndef triSurf_H
#define triSurf_H

#include "triSurfPoints.H"
#include "triSurfFacets.H"
#include "triSurfFeatureEdges.H"
#include "triSurfAddressing.H"
#include <map>
#include "DynList.H"
#include "labelLongList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Froward declarations
class triSurfModifier;

/*---------------------------------------------------------------------------*\
                        Class triSurf Declaration
\*---------------------------------------------------------------------------*/

class triSurf
:
    public triSurfPoints,
    public triSurfFacets,
    public triSurfFeatureEdges,
    public triSurfAddressing
{
    // Private member functions
        void readFromFTR(const fileName&);
        void writeToFTR(const fileName&) const;

        void readFromFMS(const fileName&);
        void writeToFMS(const fileName&) const;

        inline LongList<labelledTri>& accessToFacets();
        inline geometricSurfacePatchList& accessToPatches();

        void topologyCheck();

        //- Disallow default bitwise assignment
        void operator=(const triSurf&);

public:

    // Friend classes
        friend class triSurfModifer;

    // Constructors

        //- Default construct
        triSurf();

        //- Construct from parts
        triSurf
        (
            const LongList<labelledTri>& triangles,
            const geometricSurfacePatchList& patches,
            const edgeLongList& featureEdges_,
            const pointField& points
        );

        //- Read from file
        triSurf(const fileName& fName);

    // Destructor

        ~triSurf();


    // Member Functions
        //- read and write the surface
        void readSurface(const fileName&);
        void writeSurface(const fileName&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
