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
#include "meshSurfaceEngine.H"

//#define DEBUGTriangulation

# ifdef DEBUGTriangulation
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::calculateTrianglesAndAddressing() const
{
    if( triMeshPtr_ )
        FatalErrorIn
        (
            "void meshSurfaceOptimizer::calculateTrianglesAndAddressing() const"
        ) << "Addressing is already calculated!" << abort(FatalError);

    triMeshPtr_ = new partTriMesh(*partitionerPtr_);

    # ifdef DEBUGTriangulation
    const labelList& sPoints = triMeshPtr_->meshSurfacePointLabelInTriMesh();

    forAll(sPoints, bpI)
    {
        triSurf surf;

        const label spI = sPoints[bpI];

        partTriMeshSimplex simplex(*triMeshPtr_, spI);

        triSurfModifier sMod(surf);

        pointField& pts = sMod.pointsAccess();
        pts.setSize(simplex.pts().size());
        forAll(pts, i)
            pts[i] = simplex.pts()[i];

        LongList<labelledTri>& trias = sMod.facetsAccess();
        trias.setSize(simplex.triangles().size());
        forAll(trias, i)
        {
            const triFace& t = simplex.triangles()[i];
            trias[i] = labelledTri(t[0], t[1], t[2], 0);
        }

        sMod.patchesAccess().setSize(1);
        sMod.patchesAccess()[0].name() = "bnd";
        sMod.patchesAccess()[0].geometricType() = "patch";

        fileName sName("bndPointSimplex_");
        sName += help::scalarToText(bpI);
        sName += ".stl";
        surf.writeSurface(sName);
    }
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
