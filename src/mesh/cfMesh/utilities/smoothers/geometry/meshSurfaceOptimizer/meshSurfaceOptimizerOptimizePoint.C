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

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "meshOctree.H"
#include "triangle.H"
#include "plane.H"
#include "surfaceOptimizer.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::nodeDisplacementLaplacian
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    if( vertexType_[bpI] & LOCKED )
        return;

    const point newP = newPositionLaplacian(bpI, transformIntoPlane);

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    surfaceModifier.moveBoundaryVertex(bpI, newP);
}

void meshSurfaceOptimizer::nodeDisplacementLaplacianFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    if( vertexType_[bpI] & LOCKED )
        return;

    const point newP = newPositionLaplacianFC(bpI, transformIntoPlane);

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    surfaceModifier.moveBoundaryVertex(bpI, newP);
}

void meshSurfaceOptimizer::nodeDisplacementSurfaceOptimizer
(
    const label bpI,
    const scalar tol
)
{
    if( vertexType_[bpI] & LOCKED )
        return;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    # ifdef DEBUGSmooth
    Info << "Smoothing boundary node " << bpI << endl;
    Info << "Node label in the mesh is " << bPoints[bpI] << endl;
    Info << "Point coordinates " << points[bPoints[bpI]] << endl;
    # endif

    //- project vertices onto the plane
    const vector& pNormal = surfaceEngine_.pointNormals()[bpI];
    if( magSqr(pNormal) < VSMALL )
        return;

    const plane pl(points[bPoints[bpI]], pNormal);

    DynList<point> pts;
    DynList<triFace> trias;
    vector vecX, vecY;
    if( !transformIntoPlane(bpI, pl, vecX, vecY, pts, trias) )
    {
        Warning << "Cannot transform into plane" << endl;
        return;
    }

    surfaceOptimizer so(pts, trias);
    point newPoint = so.optimizePoint(tol);

    const point newP
    (
        points[bPoints[bpI]] +
        vecX * newPoint.x() +
        vecY * newPoint.y()
    );

    meshSurfaceEngineModifier sm(surfaceEngine_);
    sm.moveBoundaryVertex(bpI, newP);
}

void meshSurfaceOptimizer::edgeNodeDisplacement(const label bpI) const
{
    if( vertexType_[bpI] & LOCKED )
        return;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    const point pos = newEdgePositionLaplacian(bpI);
    const point newP = 0.5 * (pos + points[bPoints[bpI]]);

    # ifdef DEBUGSearch
    Info << "New position for point is " << newP << endl;
    # endif

    meshSurfaceEngineModifier(surfaceEngine_).moveBoundaryVertex(bpI, newP);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
