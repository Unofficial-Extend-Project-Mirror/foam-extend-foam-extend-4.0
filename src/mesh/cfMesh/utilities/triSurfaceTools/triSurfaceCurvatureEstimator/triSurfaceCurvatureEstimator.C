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

#include "triSurfaceCurvatureEstimator.H"
#include "demandDrivenData.H"

//#define DEBUGMorph

# ifdef DEBUGMorph
#include <sstream>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceCurvatureEstimator::triSurfaceCurvatureEstimator
(
    const triSurf& surface
)
:
    surface_(surface),
    edgePointCurvature_(),
    patchPositions_(),
    gaussianCurvature_(),
    meanCurvature_(),
    maxCurvature_(),
    minCurvature_(),
    maxCurvatureVector_(),
    minCurvatureVector_()
{
    calculateEdgeCurvature();
    calculateSurfaceCurvatures();
    //calculateGaussianCurvature();
    //calculateMeanCurvature();
    //calculateMinAndMaxCurvature();
}

triSurfaceCurvatureEstimator::~triSurfaceCurvatureEstimator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar triSurfaceCurvatureEstimator::edgePointCurvature(const label pI) const
{
    return edgePointCurvature_[pI];
}

scalar triSurfaceCurvatureEstimator::curvatureAtEdge(const label edgeI) const
{
    const edge& edg = surface_.edges()[edgeI];

    const scalar k1 = edgePointCurvature_[edg.start()];
    const scalar k2 = edgePointCurvature_[edg.end()];

    return (k1 + k2) / 2.0;
}

scalar triSurfaceCurvatureEstimator::gaussianCurvatureAtTriangle
(
    const label triI
) const
{
    const labelledTri& tri = surface_[triI];

    scalar curv(0.0);
    for(label i=0;i<3;++i)
        curv += gaussianCurvature_[tri[i]][patchPositions_(triI, i)];
    curv /= 3.0;

    return curv;
}

scalar triSurfaceCurvatureEstimator::meanCurvatureAtTriangle
(
    const label triI
) const
{
    const labelledTri& tri = surface_[triI];

    scalar curv(0.0);
    for(label i=0;i<3;++i)
        curv += meanCurvature_[tri[i]][patchPositions_(triI, i)];
    curv /= 3.0;

    return curv;
}

scalar triSurfaceCurvatureEstimator::maxCurvatureAtTriangle
(
    const label triI
) const
{
    const labelledTri& tri = surface_[triI];

    scalar curv(0.0);
    for(label i=0;i<3;++i)
        curv += maxCurvature_[tri[i]][patchPositions_(triI, i)];
    curv /= 3.0;

    return curv;
}

scalar triSurfaceCurvatureEstimator::minCurvatureAtTriangle
(
    const label triI
) const
{
    const labelledTri& tri = surface_[triI];

    scalar curv(0.0);
    for(label i=0;i<3;++i)
        curv += minCurvature_[tri[i]][patchPositions_(triI, i)];
    curv /= 3.0;

    return curv;
}

vector triSurfaceCurvatureEstimator::maxCurvatureVectorAtTriangle
(
    const label triI
) const
{
    const labelledTri& tri = surface_[triI];

    vector curv(vector::zero);
    for(label i=0;i<3;++i)
        curv += maxCurvatureVector_[tri[i]][patchPositions_(triI, i)];
    curv /= 3.0;

    return curv;
}

vector triSurfaceCurvatureEstimator::minCurvatureVectorAtTriangle
(
    const label triI
) const
{
    const labelledTri& tri = surface_[triI];

    vector curv(vector::zero);
    for(label i=0;i<3;++i)
        curv += minCurvatureVector_[tri[i]][patchPositions_(triI, i)];
    curv /= 3.0;

    return curv;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
