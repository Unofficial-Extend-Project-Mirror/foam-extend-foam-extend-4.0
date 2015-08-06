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
    triSurfaceCurvatureCalculator

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfaceCurvatureEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCurvatureEstimator_H
#define triSurfaceCurvatureEstimator_H

#include "triSurf.H"
#include "FRWGraph.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceCurvatureEstimator Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceCurvatureEstimator
{
    // Private data
        //- reference to triSurface
        const triSurf& surface_;

        //- curvature of points at feature edges
        scalarField edgePointCurvature_;

        //- curvatures of other points
        FRWGraph<label, 3> patchPositions_;
        List<DynList<scalar, 1> > gaussianCurvature_;
        List<DynList<scalar, 1> > meanCurvature_;
        List<DynList<scalar, 1> > maxCurvature_;
        List<DynList<scalar, 1> > minCurvature_;
        List<DynList<vector, 1> > maxCurvatureVector_;
        List<DynList<vector, 1> > minCurvatureVector_;

    // Private member functions
        //- calculate curvature of feature edges
        void calculateEdgeCurvature();

        //- calculate curvatures of other surface points
        void calculateSurfaceCurvatures();

        //- calculate Gaussian curvature
        void calculateGaussianCurvature();

        //- calculation of mean curvature
        void calculateMeanCurvature();

        //- calculation of min and max curvature
        void calculateMinAndMaxCurvature();

        //- Disallow default bitwise copy construct
        triSurfaceCurvatureEstimator(const triSurfaceCurvatureEstimator&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceCurvatureEstimator&);

public:

    // Constructors

        //- Construct from triSurface
        triSurfaceCurvatureEstimator(const triSurf& surface);

    // Destructor

        ~triSurfaceCurvatureEstimator();


    // Member Functions
        //- return curvature of points at feature edges
        scalar edgePointCurvature(const label) const;
        scalar curvatureAtEdge(const label) const;

        //- return average curvature on a surface triangle
        scalar gaussianCurvatureAtTriangle(const label) const;
        scalar meanCurvatureAtTriangle(const label) const;
        scalar maxCurvatureAtTriangle(const label) const;
        scalar minCurvatureAtTriangle(const label) const;

        vector maxCurvatureVectorAtTriangle(const label) const;
        vector minCurvatureVectorAtTriangle(const label) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
