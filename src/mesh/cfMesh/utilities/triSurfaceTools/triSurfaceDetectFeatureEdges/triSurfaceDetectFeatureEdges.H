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

Class
    triSurfaceDetectFeatureEdges

Description
    Divides the surface mesh into regions bounded by feature edges

SourceFiles
    triSurfaceDetectFeatureEdges.C
    triSurfaceDetectFeatureEdgesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDetectFeatureEdges_H
#define triSurfaceDetectFeatureEdges_H

#include "triSurf.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceDetectFeatureEdges Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDetectFeatureEdges
{
    // Private data
        //- reference to triSurf
        triSurf& surf_;

        //- detected feature edges
        List<direction> featureEdges_;

        //- angle tolerance
        scalar angleTolerance_;

    // Private member functions
        //- detect feature edges by check angle between the normals
        //- of the two neighbouring facets
        void detectFeatureEdgesAngleCriterion();

        //- Disallow default bitwise copy construct
        triSurfaceDetectFeatureEdges(const triSurfaceDetectFeatureEdges&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceDetectFeatureEdges&);

public:

    // Constructors

        //- Construct from triSurface and tolerance
        triSurfaceDetectFeatureEdges
        (
            triSurf& surface,
            const scalar angleDeviation = 45.0
        );

    // Destructor

        ~triSurfaceDetectFeatureEdges();


    // Member Functions

        //- store regions in a graph
        void detectFeatureEdges();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
