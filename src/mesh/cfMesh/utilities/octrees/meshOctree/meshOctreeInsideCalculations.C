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

#include "meshOctree.H"
#include "triSurf.H"
#include "helperFunctions.H"
#include "boundBox.H"

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool meshOctree::isPointInside(const point& p) const
{
    # ifdef DEBUGSearch
    Info << "Checking inside/outside for vertex " << p << endl;
    # endif
    const label cLabel = findLeafContainingVertex(p);

    if( cLabel >= 0 )
    {
        if( returnLeaf(cLabel).cubeType() & meshOctreeCubeBasic::INSIDE )
        {
            return true;
        }
    }

    return false;
}

void meshOctree::findEdgesInBox(const boundBox& bb, DynList<label>& edges) const
{
    DynList<const meshOctreeCube*, 256> neighbours;
    findLeavesContainedInBox(bb, neighbours);

    const pointField& points = surface_.points();
    const edgeLongList& sEdges = surface_.edges();
    const point c = (bb.min() + bb.max()) / 2.0;
    const scalar dSq = Foam::sqr(0.5 * (bb.max().x() - bb.min().x()));

    edges.clear();
    forAll(neighbours, neiI)
    {
        const meshOctreeCube& oc = *neighbours[neiI];

        if( !oc.hasContainedEdges() || !oc.isLeaf() )
            continue;

        const label ceI = oc.containedEdges();
        const VRWGraph& ce = oc.slotPtr()->containedEdges_;
        forAllRow(ce, ceI, i)
        {
            const edge& e = sEdges[ce(ceI, i)];
            const point p =
                help::nearestPointOnTheEdgeExact
                (
                    points[e[0]],
                    points[e[1]],
                    c
                );

            if( magSqr(p - c) < dSq )
                edges.append(ce(ceI, i));
        }
    }
}

void meshOctree::findTrianglesInBox
(
    const boundBox& bb,
    DynList<label>& triaLabels
) const
{
    DynList<const meshOctreeCube*, 256> neighbours;
    findLeavesContainedInBox(bb, neighbours);

    const point c = (bb.min() + bb.max()) / 2.0;
    const scalar dSq = Foam::sqr(0.5 * (bb.max().x() - bb.min().x()));

    triaLabels.clear();
    forAll(neighbours, neiI)
    {
        const meshOctreeCube& oc = *neighbours[neiI];

        if( !oc.hasContainedElements() || !oc.isLeaf() )
            continue;

        const label elI = oc.containedElements();
        const VRWGraph& ct = oc.slotPtr()->containedTriangles_;
        forAllRow(ct, elI, i)
        {
            const label tI = ct(elI, i);

            const point p = help::nearestPointOnTheTriangle(tI, surface_, c);
            if( Foam::magSqr(p - c) < dSq )
                triaLabels.append(tI);
        }
    }
}

void meshOctree::findLeavesContainedInBox
(
    const boundBox& bb,
    DynList<const meshOctreeCube*, 256>& containedCubes
) const
{
    containedCubes.clear();

    initialCubePtr_->leavesInBox(rootBox_, bb, containedCubes);
}

void meshOctree::findLeavesContainedInBox
(
    const boundBox& bb,
    labelList& containedCubes
) const
{
    DynList<const meshOctreeCube*, 256> cb;
    findLeavesContainedInBox(bb, cb);

    containedCubes.setSize(cb.size());
    label i(0);
    forAll(cb, cI)
    {
        if( cb[cI]->isLeaf() )
        {
            containedCubes[i] = cb[cI]->cubeLabel();
            ++i;
        }
    }

    containedCubes.setSize(i);
}

void meshOctree::findLeavesContainedInBox
(
    const boundBox& bb,
    DynList<label>& containedCubes
) const
{
    DynList<const meshOctreeCube*, 256> cb;
    findLeavesContainedInBox(bb, cb);

    containedCubes.clear();

    forAll(cb, cI)
    {
        if( cb[cI]->isLeaf() )
            containedCubes.append(cb[cI]->cubeLabel());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
