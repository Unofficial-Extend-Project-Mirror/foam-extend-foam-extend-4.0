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

#include "triSurfaceClassifyEdges.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "triSurf.H"
#include "meshOctree.H"
#include "labelPair.H"

#ifdef USE_OMP
#include <omp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceClassifyEdges::checkOrientation()
{
    const triSurf& surf = octree_.surface();
    const boundBox& rootBox = octree_.rootBox();
    const pointField& points = surf.points();
    const VRWGraph& facetEdges = surf.facetEdges();
    const VRWGraph& edgeFacets = surf.edgeFacets();

    facetOrientation_.setSize(surf.size());

    //- sort all surface facets into groups consisting of facets with consistent
    //- orientation. Do not cross non-manifold edges
    labelLongList orientationGroup(surf.size(), -1);
    label nGroups(0);

    forAll(surf, triI)
    {
        if( orientationGroup[triI] != -1 )
            continue;

        orientationGroup[triI] = nGroups;
        labelLongList front;
        front.append(triI);

        while( front.size() != 0 )
        {
            const label tLabel = front.removeLastElement();

            const labelledTri& facet = surf[tLabel];

            forAll(facet, eI)
            {
                const label edgeI = facetEdges(tLabel, eI);

                if( edgeFacets.sizeOfRow(edgeI) != 2 )
                    continue;

                forAllRow(edgeFacets, edgeI, efI)
                {
                    const label neiFacetI = edgeFacets(edgeI, efI);

                    if( orientationGroup[neiFacetI] != -1 )
                        continue;
                    if( neiFacetI == tLabel )
                        continue;

                    const labelledTri& neiFacet = surf[neiFacetI];

                    //- check the orientation of triangles at this edge
                    //- check the sign of the angle
                    //- if the orientation  is not consistent
                    DynList<labelPair, 2> sharedIndices;
                    forAll(facet, i)
                    {
                        forAll(neiFacet, j)
                        {
                            if( facet[i] == neiFacet[j] )
                                sharedIndices.append(labelPair(i, j));
                        }
                    }

                    if( sharedIndices.size() == 2 )
                    {
                        const labelPair& pair0 = sharedIndices[0];
                        const labelPair& pair1 = sharedIndices[1];
                        if( ((pair0.first() + 1) % 3) == pair1.first() )
                        {
                            if( (pair1.second() + 1) % 3 == pair0.second() )
                            {
                                orientationGroup[neiFacetI] = nGroups;
                                front.append(neiFacetI);
                            }
                        }
                        else
                        {
                            if( (pair0.second() + 1) % 3 == pair1.second() )
                            {
                                orientationGroup[neiFacetI] = nGroups;
                                front.append(neiFacetI);
                            }
                        }
                    }
                }
            }
        }

        ++nGroups;
    }

    Info << "Found " << nGroups
         << " groups of triangles with consistent orientation" << endl;

    //- find the octree leaves containing each triangle
    VRWGraph triangleInLeaves(surf.size());
    labelLongList ntl(surf.size(), 0);

    DynList<label> helper;
    for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
    {
        helper.clear();
        octree_.containedTriangles(leafI, helper);

        forAll(helper, i)
            ++ntl[helper[i]];
    }

    forAll(ntl, triI)
        triangleInLeaves.setRowSize(triI, ntl[triI]);

    ntl = 0;
    for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
    {
        helper.clear();
        octree_.containedTriangles(leafI, helper);

        forAll(helper, i)
        {
            const label triI = helper[i];

            triangleInLeaves(triI, ntl[triI]++) = leafI;
        }
    }

    //- check the orientation of all facets in a group and collect their votes
    DynList<labelPair> groupVotes;
    groupVotes.setSize(nGroups);
    groupVotes = labelPair(0, 0);

    # ifdef USE_OMP
    # pragma omp parallel if( surf.size() > 1000 ) private(helper)
    # endif
    {
        DynList<labelPair> localVotes;
        localVotes.setSize(nGroups);
        localVotes = labelPair(0, 0);

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(orientationGroup, triI)
        {
            const labelledTri& tri = surf[triI];
            const point c = tri.centre(points);
            vector n = tri.normal(points);
            const scalar magN = mag(n);

            if( magN < VSMALL )
                continue;

            n /= magN;

            //- find the OUTSIDE octree cubes in the vicinity of the triangle
            //- and check the orientation of the triangle
            forAllRow(triangleInLeaves, triI, tlI)
            {
                const label leafI = triangleInLeaves(triI, tlI);

                octree_.findAllLeafNeighbours(leafI, helper);

                forAll(helper, i)
                {
                    const label leafJ = helper[i];

                    const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafJ);

                    if( oc.cubeType() & meshOctreeCubeBasic::OUTSIDE )
                    {
                        const scalar length = 3.0 * oc.size(rootBox);

                        point pMin, pMax;
                        oc.cubeBox(rootBox, pMin, pMax);

                        const boundBox bb(pMin, pMax);

                        //- check whether the ray casted from
                        //- the centre of the triangle intersects the cube
                        const point endPos = c + length * n;
                        const point endNeg = c - length * n;

                        if( help::boundBoxLineIntersection(c, endPos, bb) )
                        {
                            //- found an intersection in the positive direction
                            ++localVotes[orientationGroup[triI]].first();
                        }
                        else if( help::boundBoxLineIntersection(c, endNeg, bb) )
                        {
                            //- found an intersection in the negative direction
                            ++localVotes[orientationGroup[triI]].second();
                        }
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp critical(grouping)
        # endif
        {
            forAll(localVotes, groupI)
            {
                groupVotes[groupI].first() += localVotes[groupI].first();
                groupVotes[groupI].second() += localVotes[groupI].second();
            }
        }
    }

    Info << "Before determining of orientation" << endl;

    //- determine whether a group is oriented outward or inward
    List<direction> outwardGroup(nGroups, direction(0));

    forAll(groupVotes, groupI)
    {
        if( groupVotes[groupI].first() > groupVotes[groupI].second() )
        {
            outwardGroup[groupI] = 1;
        }
        else if( groupVotes[groupI].first() < groupVotes[groupI].second() )
        {
            outwardGroup[groupI] = 2;
        }
    }

    Info << "Here setting orientation" << endl;
    //- Finally, set the orientation of the normal
    facetOrientation_.setSize(surf.size());
    forAll(facetOrientation_, triI)
        facetOrientation_[triI] = outwardGroup[orientationGroup[triI]];
}

void triSurfaceClassifyEdges::classifyEdgesTypes()
{
    const triSurf& surf = octree_.surface();
    const pointField& points = surf.points();
    const VRWGraph& edgeFacets = surf.edgeFacets();
    const edgeLongList& edges = surf.edges();
    const VRWGraph& pointEdges = surf.pointEdges();
    const edgeLongList& featureEdges = surf.featureEdges();

    edgeTypes_.setSize(edgeFacets.size());
    edgeTypes_ = NONE;

    //- set feature edge types
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(featureEdges, feI)
    {
        const edge& e = featureEdges[feI];

        forAllRow(pointEdges, e.start(), peI)
        {
            if( edges[pointEdges(e.start(), peI)] == e )
                edgeTypes_[pointEdges(e.start(), peI)] |= FEATUREEDGE;
        }
    }

    //- Finally, check the edges and assign flags
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(edgeFacets, edgeI)
    {
        if( edgeFacets.sizeOfRow(edgeI) != 2 )
        {
            //- this is not a manifold edge
            edgeTypes_[edgeI] = NONE;
            continue;
        }

        //- surface is a manifold at this edge
        const label f0 = edgeFacets(edgeI, 0);
        const label f1 = edgeFacets(edgeI, 1);

        const labelledTri& tri0 = surf[f0];
        const labelledTri& tri1 = surf[f1];

        if( tri0.region() != tri1.region() )
            edgeTypes_[edgeI] |= FEATUREEDGE;

        label apexNode(-1);
        forAll(tri1, pI)
        {
            bool found(false);
            forAll(tri0, pJ)
            {
                if( tri0[pJ] == tri1[pI] )
                {
                    found = true;
                    break;
                }

                if( found )
                    break;
            }

            if( !found )
            {
                apexNode = tri1[pI];
                break;
            }
        }

        const tetrahedron<point, point> tet
        (
            points[tri0[0]],
            points[tri0[1]],
            points[tri0[2]],
            points[apexNode]
        );

        const scalar vol = tet.mag();

        if( facetOrientation_[f0] == 1 )
        {
            //- facet is outward oriented
            if( vol < -VSMALL )
            {
                //- this is a convex edge
                edgeTypes_[edgeI] |= CONVEXEDGE;
            }
            else if( vol > VSMALL )
            {
                //- this is a concave edge
                edgeTypes_[edgeI] |= CONCAVEEDGE;
            }
            else
            {
                edgeTypes_[edgeI] |= FLATSURFACEEDGE;
            }
        }
        else if( facetOrientation_[f0] == 2 )
        {
            //- facet is inward oriented
            if( vol > VSMALL )
            {
                //- this is a convex edge
                edgeTypes_[edgeI] |= CONVEXEDGE;
            }
            else if( vol < -VSMALL )
            {
                //- this is a concave edge
                edgeTypes_[edgeI] |= CONCAVEEDGE;
            }
            else
            {
                edgeTypes_[edgeI] |= FLATSURFACEEDGE;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
