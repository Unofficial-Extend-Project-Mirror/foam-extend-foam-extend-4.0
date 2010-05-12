/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "DynamicList.H"
#include "demandDrivenData.H"
#include "SortableList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Returns edgeI between two points.
Foam::label primitiveMesh::getEdge
(
    List<DynamicList<label> >& pe,
    DynamicList<edge>& es,

    const label pointI,
    const label nextPointI
)
{
    // Find connection between pointI and nextPointI
    forAll(pe[pointI], ppI)
    {
        label eI = pe[pointI][ppI];

        const edge& e = es[eI];

        if (e.start() == nextPointI || e.end() == nextPointI)
        {
            return eI;
        }
    }

    // Make new edge.
    label edgeI = es.size();
    pe[pointI].append(edgeI);
    pe[nextPointI].append(edgeI);
    if (pointI < nextPointI)
    {
        es.append(edge(pointI, nextPointI));
    }
    else
    {
        es.append(edge(nextPointI, pointI));
    }
    return edgeI;
}


void primitiveMesh::calcEdges(const bool doFaceEdges) const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcEdges(const bool) : "
            << "calculating edges, pointEdges and optionally faceEdges"
            << endl;
    }

    // It is an error to attempt to recalculate edges
    // if the pointer is already set
    if ((edgesPtr_ || pePtr_) || (doFaceEdges && fePtr_))
    {
        FatalErrorIn("primitiveMesh::calcEdges(const bool) const")
            << "edges or pointEdges or faceEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        // ALGORITHM:
        // Go through the faces list. Search pointEdges for existing edge.
        // If not found create edge and add to pointEdges.
        // In second loop sort edges in order of increasing neighbouring
        // point.
        // This algorithm replaces the one using pointFaces which used more
        // allocations but less memory and was on practical cases
        // quite a bit slower.

        const faceList& fcs = faces();

        // Size up lists
        // ~~~~~~~~~~~~~

        // Estimate pointEdges storage
        List<DynamicList<label> > pe(nPoints());
        forAll(pe, pointI)
        {
            pe[pointI].setSize(primitiveMesh::edgesPerPoint_);
        }

        // Estimate edges storage
        DynamicList<edge> es(pe.size()*primitiveMesh::edgesPerPoint_/2);

        // Estimate faceEdges storage
        if (doFaceEdges)
        {
            fePtr_ = new labelListList(fcs.size());
            labelListList& faceEdges = *fePtr_;
            forAll(fcs, faceI)
            {
                faceEdges[faceI].setSize(fcs[faceI].size());
            }
        }


        // Find consecutive face points in edge list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Edges on boundary faces
        label nExtEdges = 0;
        // Edges using no boundary point
        nInternal0Edges_ = 0;
        // Edges using 1 boundary point
        label nInt1Edges = 0;
        // Edges using two boundary points but not on boundary face:
        // edges.size()-nExtEdges-nInternal0Edges_-nInt1Edges

        // Ordering is different if points are ordered.
        if (nInternalPoints_ == -1)
        {
            // No ordering. No distinction between types.
            forAll(fcs, faceI)
            {
                const face& f = fcs[faceI];

                forAll(f, fp)
                {
                    label pointI = f[fp];
                    label nextPointI = f[f.fcIndex(fp)];

                    label edgeI = getEdge(pe, es, pointI, nextPointI);

                    if (doFaceEdges)
                    {
                        (*fePtr_)[faceI][fp] = edgeI;
                    }
                }
            }
            // Assume all edges internal
            nExtEdges = 0;
            nInternal0Edges_ = es.size();
            nInt1Edges = 0;
        }
        else
        {
            // 1. Do external faces first. This creates external edges.
            for (label faceI = nInternalFaces_; faceI < fcs.size(); faceI++)
            {
                const face& f = fcs[faceI];

                forAll(f, fp)
                {
                    label pointI = f[fp];
                    label nextPointI = f[f.fcIndex(fp)];

                    label oldNEdges = es.size();
                    label edgeI = getEdge(pe, es, pointI, nextPointI);

                    if (es.size() > oldNEdges)
                    {
                        nExtEdges++;
                    }
                    if (doFaceEdges)
                    {
                        (*fePtr_)[faceI][fp] = edgeI;
                    }
                }
            }

            // 2. Do internal faces. This creates internal edges.
            for (label faceI = 0; faceI < nInternalFaces_; faceI++)
            {
                const face& f = fcs[faceI];

                forAll(f, fp)
                {
                    label pointI = f[fp];
                    label nextPointI = f[f.fcIndex(fp)];

                    label oldNEdges = es.size();
                    label edgeI = getEdge(pe, es, pointI, nextPointI);

                    if (es.size() > oldNEdges)
                    {
                        if (pointI < nInternalPoints_)
                        {
                            if (nextPointI < nInternalPoints_)
                            {
                                nInternal0Edges_++;
                            }
                            else
                            {
                                nInt1Edges++;
                            }
                        }
                        else
                        {
                            if (nextPointI < nInternalPoints_)
                            {
                                nInt1Edges++;
                            }
                            else
                            {
                                // Internal edge with two points on boundary
                            }
                        }
                    }
                    if (doFaceEdges)
                    {
                        (*fePtr_)[faceI][fp] = edgeI;
                    }
                }
            }
        }


        // For unsorted meshes the edges will be in order of occurrence inside
        // the faces. For sorted meshes the first nExtEdges will be the external
        // edges.

        if (nInternalPoints_ != -1)
        {
            nInternalEdges_ = es.size()-nExtEdges;
            nInternal1Edges_ = nInternal0Edges_+nInt1Edges;

            //Pout<< "Edge overview:" << nl
            //    << "    total number of edges           : " << es.size()
            //    << nl
            //    << "    boundary edges                  : " << nExtEdges
            //    << nl
            //    << "    edges using no boundary point   : "
            //    << nInternal0Edges_
            //    << nl
            //    << "    edges using one boundary points : " << nInt1Edges
            //   << nl
            //    << "    edges using two boundary points : "
            //    << es.size()-nExtEdges-nInternal0Edges_-nInt1Edges << endl;
        }


        // Like faces sort edges in order of increasing neigbouring point.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Automatically if points are sorted into internal and external points
        // the edges will be sorted into
        // - internal point to internal point
        // - internal point to external point
        // - external point to external point
        // The problem is that the last one mixes external edges with internal
        // edges with two points on the boundary.


        // Map to sort into new upper-triangular order
        labelList oldToNew(es.size());
        if (debug)
        {
            oldToNew = -1;
        }

        // start of edges with 0 boundary points
        label internal0EdgeI = 0;

        // start of edges with 1 boundary points
        label internal1EdgeI = nInternal0Edges_;

        // start of edges with 2 boundary points
        label internal2EdgeI = nInternal1Edges_;

        // start of external edges
        label externalEdgeI = nInternalEdges_;


        // To sort neighbouring points in increasing order. Defined outside
        // for optimisation reasons: if all connectivity size same will need
        // no reallocations
        SortableList<label> nbrPoints(primitiveMesh::edgesPerPoint_);

        forAll(pe, pointI)
        {
            const DynamicList<label>& pEdges = pe[pointI];

            nbrPoints.setSize(pEdges.size());

            forAll(pEdges, i)
            {
                const edge& e = es[pEdges[i]];

                label nbrPointI = e.otherVertex(pointI);

                if (nbrPointI < pointI)
                {
                    nbrPoints[i] = -1;
                }
                else
                {
                    nbrPoints[i] = nbrPointI;
                }
            }
            nbrPoints.sort();


            if (nInternalPoints_ == -1)
            {
                // Sort all upper-triangular
                forAll(nbrPoints, i)
                {
                    if (nbrPoints[i] != -1)
                    {
                        label edgeI = pEdges[nbrPoints.indices()[i]];

                        oldToNew[edgeI] = internal0EdgeI++;
                    }
                }        
            }
            else
            {
                if (pointI < nInternalPoints_)
                {
                    forAll(nbrPoints, i)
                    {
                        label nbrPointI = nbrPoints[i];

                        label edgeI = pEdges[nbrPoints.indices()[i]];

                        if (nbrPointI != -1)
                        {
                            if (edgeI < nExtEdges)
                            {
                                // External edge
                                oldToNew[edgeI] = externalEdgeI++;
                            }
                            else if (nbrPointI < nInternalPoints_)
                            {
                                // Both points inside
                                oldToNew[edgeI] = internal0EdgeI++;
                            }
                            else
                            {
                                // One points inside, one outside
                                oldToNew[edgeI] = internal1EdgeI++;
                            }
                        }
                    }
                }
                else
                {
                    forAll(nbrPoints, i)
                    {
                        label nbrPointI = nbrPoints[i];

                        label edgeI = pEdges[nbrPoints.indices()[i]];

                        if (nbrPointI != -1)
                        {
                            if (edgeI < nExtEdges)
                            {
                                // External edge
                                oldToNew[edgeI] = externalEdgeI++;
                            }
                            else if (nbrPointI < nInternalPoints_)
                            {
                                // Not possible!
                                FatalErrorIn("primitiveMesh::calcEdges(..)")
                                    << abort(FatalError);
                            }
                            else
                            {
                                // Both points outside
                                oldToNew[edgeI] = internal2EdgeI++;
                            }
                        }
                    }
                }
            }
        }


        if (debug)
        {
            label edgeI = findIndex(oldToNew, -1);

            if (edgeI != -1)
            {
                const edge& e = es[edgeI];

                FatalErrorIn("primitiveMesh::calcEdges(..)")
                    << "Did not sort edge " << edgeI << " points:" << e
                    << " coords:" << points()[e[0]] << points()[e[1]]
                    << endl
                    << "Current buckets:" << endl
                    << "    internal0EdgeI:" << internal0EdgeI << endl
                    << "    internal1EdgeI:" << internal1EdgeI << endl
                    << "    internal2EdgeI:" << internal2EdgeI << endl
                    << "    externalEdgeI:" << externalEdgeI << endl
                    << abort(FatalError);
            }
        }



        // Renumber and transfer.

        // Edges
        edgesPtr_ = new edgeList(es.size());
        edgeList& edges = *edgesPtr_;
        forAll(es, edgeI)
        {
            edges[oldToNew[edgeI]] = es[edgeI];
        }

        // pointEdges
        pePtr_ = new labelListList(nPoints());
        labelListList& pointEdges = *pePtr_;
        forAll(pe, pointI)
        {
            DynamicList<label>& pEdges = pe[pointI];
            inplaceRenumber(oldToNew, pEdges);
            pEdges.shrink();
            pointEdges[pointI].transfer(pEdges);
            Foam::sort(pointEdges[pointI]);
            pEdges.clear();
        }

        // faceEdges
        if (doFaceEdges)
        {
            labelListList& faceEdges = *fePtr_;
            forAll(faceEdges, faceI)
            {
                inplaceRenumber(oldToNew, faceEdges[faceI]);
            }
        }
    }
}


void primitiveMesh::calcOrderedEdges() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcOrderedEdges() : "
            << "calculating edges and faceEdges"
            << endl;
    }

    // It is an error to attempt to recalculate ordered edges
    // if the pointer is already set
    if (orderedEdgesPtr_ )
    {
        FatalErrorIn("primitiveMesh::calcOrderedEdges() const")
            << "edges and faceEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        // ALGORITHM:
        // Go through the pointFace list.  Go through the list of faces for that
        // point and ask for edges.  If the edge has got the point in question
        // AND the second point in the edge is larger than the first, add the
        // edge to the list.  At the same time, add the edge label to the list
        // of edges for the current face (faceEdges) and log the other face as
        // the neighbour of this face.

        const faceList& f = faces();

        const labelListList& pf = pointFaces();

        labelListList fe(nFaces());

        // count the maximum number of edges
        label maxEdges = 0;

        // create a list of edges for each face and store for efficiency
        edgeListList edgesOfFace(nFaces());

        forAll (f, faceI)
        {
            edgesOfFace[faceI] = f[faceI].edges();

            maxEdges += f[faceI].nEdges();

            labelList& curFE = fe[faceI];

            curFE.setSize(f[faceI].nEdges());

            forAll (curFE, curFEI)
            {
                curFE[curFEI] = -1;
            }
        }

        // EDGE CALCULATION

        orderedEdgesPtr_ = new edgeList(maxEdges);
        edgeList& e = *orderedEdgesPtr_;
        label nEdges = 0;

        forAll (pf, pointI)
        {
            const labelList& curFaces = pf[pointI];

            // create a list of labels to keep the neighbours that
            // have already been added
            DynamicList<label, edgesPerPoint_> addedNeighbours;
            DynamicList<DynamicList<label, edgesPerPoint_> >
                faceGivingNeighbour;
            DynamicList<DynamicList<label, edgesPerPoint_> >
                edgeOfFaceGivingNeighbour;

            forAll (curFaces, faceI)
            {
                // get the edges
                const edgeList& fEdges = edgesOfFace[curFaces[faceI]];

                // for every edge
                forAll(fEdges, edgeI)
                {
                    const edge& ends = fEdges[edgeI];

                    // does the edge has got the point in question
                    bool found = false;
                    label secondPoint = -1;

                    if (ends.start() == pointI)
                    {
                        found = true;
                        secondPoint = ends.end();
                    }

                    if (ends.end() == pointI)
                    {
                        found = true;
                        secondPoint = ends.start();
                    }

                    // if the edge has got the point and second label is larger
                    // than first, it is a candidate for adding
                    if (found && (secondPoint > pointI))
                    {
                        // check if the edge has already been added
                        bool added = false;

                        forAll (addedNeighbours, eopI)
                        {
                            if (secondPoint == addedNeighbours[eopI])
                            {
                                // Edge is already added. New face sharing it
                                added = true;

                                // Remember the face and edge giving neighbour
                                faceGivingNeighbour[eopI].append
                                    (curFaces[faceI]);

                                edgeOfFaceGivingNeighbour[eopI].append(edgeI);

                                break;
                            }
                        }

                        // If not added, add the edge to the list
                        if (!added)
                        {
                            addedNeighbours.append(secondPoint);

                            // Remember the face and subShape giving neighbour
                            faceGivingNeighbour(addedNeighbours.size() - 1)
                                .append(curFaces[faceI]);
                            edgeOfFaceGivingNeighbour
                                (addedNeighbours.size() - 1).append(edgeI);
                        }
                    }
                }
            }

            // All edges for the current point found. Before adding them to the
            // list, it is necessary to sort them in the increasing order of the
            // neighbouring point.

            // Make real list out of SLList to simplify the manipulation.
            // Also, make another list to "remember" how the original list was
            // reshuffled.
            labelList shuffleList(addedNeighbours.size());

            forAll (shuffleList, i)
            {
                shuffleList[i] = i;
            }

            // Use a simple sort to sort the addedNeighbours list.
            //  Other two lists mimic the same behaviour
            label i, j, a, b;

            for (j = 1; j <= addedNeighbours.size() - 1; j++)
            {
                a = addedNeighbours[j];
                b = shuffleList[j];

                i = j - 1;

                while (i >= 0 && addedNeighbours[i] > a)
                {
                    addedNeighbours[i + 1] = addedNeighbours[i];
                    shuffleList[i + 1] = shuffleList[i];
                    i--;
                }

                addedNeighbours[i + 1] = a;
                shuffleList[i + 1] = b;
            }

            labelList reshuffleList(shuffleList.size());

            forAll(shuffleList, i)
            {
                reshuffleList[shuffleList[i]] = i;
            }

            // Reshuffle other lists

            labelListList fgn(faceGivingNeighbour.size());

            forAll (faceGivingNeighbour, i)
            {
                fgn[reshuffleList[i]].transfer(faceGivingNeighbour[i].shrink());
            }

            labelListList eofgn(edgeOfFaceGivingNeighbour.size());

            forAll (edgeOfFaceGivingNeighbour, i)
            {
                eofgn[reshuffleList[i]].transfer
                (
                    edgeOfFaceGivingNeighbour[i].shrink()
                );
            }

            // adding the edges
            forAll(addedNeighbours, edgeI)
            {
                const labelList& curFgn = fgn[edgeI];
                const labelList& curEofgn = eofgn[edgeI];

                forAll (curFgn, fgnI)
                {
                    fe[curFgn[fgnI]][curEofgn[fgnI]] = nEdges;
                }

                e[nEdges] = edge(pointI, addedNeighbours[edgeI]);
                nEdges++;
            }
        }

        // reset the size
        e.setSize(nEdges);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const edgeList& primitiveMesh::edges() const
{
    if (!edgesPtr_)
    {
        calcEdges(true);
    }

    return *edgesPtr_;
}


const labelListList& primitiveMesh::pointEdges() const
{
    if (!pePtr_)
    {
        //// Invert edges
        //pePtr_ = new labelListList(nPoints());
        //invertManyToMany(nPoints(), edges(), *pePtr_);
        calcEdges(true);
    }

    return *pePtr_;
}


const labelListList& primitiveMesh::faceEdges() const
{
    if (!fePtr_)
    {
        calcEdges(true);
    }

    return *fePtr_;
}


// Ordered edge addressing.  HJ, 6/Jan/2009
const edgeList& primitiveMesh::orderedEdges() const
{
    if (!orderedEdgesPtr_)
    {
        calcOrderedEdges();
    }

    return *orderedEdgesPtr_;
}


void primitiveMesh::clearOutEdges()
{
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(fePtr_);

    deleteDemandDrivenData(orderedEdgesPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
