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

#include "voronoiMeshExtractor.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGVoronoi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::createAddressing() const
{
    if( pointEdgesPtr_ || edgeTetsPtr_ || boundaryEdgePtr_ || edgesPtr_ )
        return;

    pointEdgesPtr_ = new VRWGraph(tetCreator_.tetPoints().size());
    VRWGraph& pointEdges = *pointEdgesPtr_;

    edgeTetsPtr_ = new VRWGraph();
    VRWGraph& edgeTets = *edgeTetsPtr_;

    boundaryEdgePtr_ = new boolList();
    boolList& boundaryEdges = *boundaryEdgePtr_;

    edgesPtr_ = new LongList<edge>();
    LongList<edge>& edges = *edgesPtr_;

    //- create edges and edgeTets
    const LongList<partTet>& tets = tetCreator_.tets();

    VRWGraph pointTets;
    pointTets.reverseAddressing(tets);

    forAll(tets, tetI)
    {
        const FixedList<edge, 6> tetEdges = tets[tetI].edges();

        forAll(tetEdges, eI)
        {
            const edge& e = tetEdges[eI];
            const label start = e.start();

            const row endTets = pointTets[e.end()];

            bool store(true);

            DynList<label> eTets;
            forAllRow(pointTets, start, ptI)
            {
                const label tetJ = pointTets(start, ptI);

                if( !endTets.contains(tetJ) )
                    continue;

                if( tetJ < tetI )
                {
                    store = false;
                    break;
                }

                eTets.append(tetJ);
            }

            if( store )
            {
                edgeTets.appendList(eTets);

                edges.append(e);
            }
        }
    }

    # ifdef DEBUGVoronoi
    Info << "Edge tets " << edgeTets << endl;
    # endif

    //- calculate point-edges addressing
    pointEdges.reverseAddressing(edges);

    //- sort edge-tets in circular order
    boundaryEdges.setSize(edgeTets.size());
    boundaryEdges = false;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeTets, edgeI)
    {
        const edge& e = edges[edgeI];
        row eTets = edgeTets[edgeI];

        # ifdef DEBUGVoronoi
        Info << "Edge " << edgeI << " has points " << e << endl;
        # endif

        forAll(eTets, tetI)
        {
            const partTet& pt = tets[eTets[tetI]];

            # ifdef DEBUGVoronoi
            Info << "Checking tet " << eTets[tetI] << " points " << pt << endl;
            # endif

            //- find the face shared with the neighbour
            const FixedList<edge, 6> tetEdges = pt.edges();

            label searchPoint(-1);
            forAll(tetEdges, eI)
            {
                if( tetEdges[eI] == e )
                {
                    if( tetEdges[eI].start() == e.start() )
                    {
                        searchPoint = pt[sameOrientation_[eI]];
                    }
                    else
                    {
                        searchPoint = pt[oppositeOrientation_[eI]];
                    }

                    break;
                }
            }

            if( searchPoint < 0 )
                FatalErrorIn
                (
                    "void voronoiMeshExtractor::createAddressing() const"
                ) << " invalid search point " << abort(FatalError);

            bool found(false);
            forAll(eTets, i)
            {
                if( tetI == i )
                    continue;

                const partTet& ptNei = tets[eTets[i]];

                const label pos = ptNei.whichPosition(searchPoint);

                if( pos < 0 )
                    continue;

                if( tetI < eTets.size()-1 )
                {
                    const label add  = eTets[tetI+1];
                    eTets[tetI+1] = eTets[i];
                    eTets[i] = add;
                }
                found = true;
                break;
            }

            if( !found )
                boundaryEdges[edgeI] = true;
        }
    }
}

const VRWGraph& voronoiMeshExtractor::pointEdges() const
{
    if( !pointEdgesPtr_ )
        createAddressing();

    return *pointEdgesPtr_;
}

const LongList<edge>& voronoiMeshExtractor::edges() const
{
    if( !edgesPtr_ )
        createAddressing();

    return *edgesPtr_;
}

const VRWGraph& voronoiMeshExtractor::edgeTets() const
{
    if( !edgeTetsPtr_ )
        createAddressing();

    return *edgeTetsPtr_;
}

const boolList& voronoiMeshExtractor::boundaryEdge() const
{
    if( !boundaryEdgePtr_ )
        createAddressing();

    return *boundaryEdgePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
