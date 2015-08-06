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

#include "triSurfaceExtrude2DEdges.H"
#include "triSurfModifier.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceExtrude2DEdges::triSurfaceExtrude2DEdges(const triSurf& surface)
:
    surf_(surface)
{}

triSurfaceExtrude2DEdges::~triSurfaceExtrude2DEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceExtrude2DEdges::extrudeSurface(triSurf& newSurf) const
{
    triSurfModifier sMod(newSurf);

    //- set patches
    geometricSurfacePatchList& patches = sMod.patchesAccess();

    patches.setSize(1);
    patches[0].name() = "patch0";
    patches[0].geometricType() = "patch";

    //- check if the edges are in the x-y plane
    const pointField& sPoints = surf_.points();
    const boundBox bb(sPoints);

    if( Foam::mag(bb.max().z() - bb.min().z()) > SMALL )
        FatalErrorIn
        (
            "void triSurfaceExtrude2DEdges::extrudeSurface(triSurf&) const"
        ) << "Cannot extrude edges which are not in the x-y plane!"
          << exit(FatalError);

    //- copy points
    pointField& pts = sMod.pointsAccess();
    pts.setSize(2 * sPoints.size());

    const label nOffset = sPoints.size();
    const scalar zOffset = 0.1 * bb.mag();

    forAll(sPoints, pI)
    {
        pts[pI] = pts[pI+nOffset] = sPoints[pI];
        pts[pI+sPoints.size()].z() += zOffset;
    }

    //- create triangles from feature edges
    LongList<labelledTri>& triangles = sMod.facetsAccess();
    const edgeLongList& edges = surf_.featureEdges();

    triangles.setSize(2 * edges.size());
    forAll(edges, eI)
    {
        const edge& e = edges[eI];
        const label tI = 2 * eI;
        triangles[tI] = labelledTri(e[0], e[1], e[1]+nOffset, 0);
        triangles[tI + 1] = labelledTri(e[0], e[1]+nOffset, e[0]+nOffset, 0);
    }
}

const triSurf* triSurfaceExtrude2DEdges::extrudeSurface() const
{
    triSurf* sPtr = new triSurf();

    extrudeSurface(*sPtr);

    return sPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
