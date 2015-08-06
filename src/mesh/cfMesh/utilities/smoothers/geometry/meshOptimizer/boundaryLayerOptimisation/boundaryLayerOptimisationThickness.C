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
#include "boundaryLayerOptimisation.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "labelledScalar.H"
#include "polyMeshGenAddressing.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::hairEdgesAtBndFace
(
    const label cellI,
    const label baseFaceI,
    DynList<edge>& hairEdges
) const
{
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

        //- check cell topology
    DynList<edge, 48> edges;
    DynList<DynList<label, 2>, 48> edgeFaces;
    DynList<DynList<label, 10>, 24> faceEdges;
    faceEdges.setSize(c.size());
    label baseFace(-1);
    forAll(c, fI)
    {
        if( c[fI] == baseFaceI )
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = edges.containsAtPosition(e);

            if( pos < 0 )
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos+1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }

    const face& bf = faces[c[baseFace]];
    hairEdges.setSize(bf.size());

    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][bf.rcIndex(pI)];

        if( edgeFaces[nextEdge].size() != 2 || edgeFaces[prevEdge].size() != 2 )
            break;

        //- find the face attached to the edge after the current point
        label otherNextFace = edgeFaces[nextEdge][0];
        if( otherNextFace == baseFace )
            otherNextFace = edgeFaces[nextEdge][1];

        //- find the face attached to the edge before the current point
        label otherPrevFace = edgeFaces[prevEdge][0];
        if( otherPrevFace == baseFace )
            otherPrevFace = edgeFaces[prevEdge][1];

        label commonEdge;
        for(commonEdge=0;commonEdge<edges.size();++commonEdge)
            if(
                edgeFaces[commonEdge].contains(otherNextFace) &&
                edgeFaces[commonEdge].contains(otherPrevFace)
            )
                break;

        if( commonEdge == edges.size() )
            break;

        //- there exists a common edge which shall be used as a hair
        if( edges[commonEdge].start() == bf[pI] )
        {
            hairEdges[pI] = edges[commonEdge];
        }
        else
        {
            hairEdges[pI] = edges[commonEdge].reverseEdge();
        }
    }
}

scalar boundaryLayerOptimisation::calculateThickness
(
    const label heI,
    const label heJ
) const
{
    const pointFieldPMG& points = mesh_.points();

    //- references to hair edges
    const edge& he = hairEdges_[heI];
    const edge& nhe = hairEdges_[heJ];

    //- distance vector between the surface points of hair edges
    const point& sp = points[he[0]];
    const point& ep = points[nhe[0]];
    const vector dv = ep - sp;
    const scalar magDv = mag(dv);

    //- calculate layer thickness
    const scalar currThickness = he.mag(points);
    scalar retThickness = currThickness;

    const scalar currNeiThickness = nhe.mag(points);
    scalar suggestedNeiThickness = currNeiThickness;

    //- calculate layer height at the current point
    const point npAlpha = help::nearestPointOnTheEdge(sp, ep, points[he[1]]);
    const scalar currHeight = mag(npAlpha - points[he[1]]);
    scalar retHeight = currHeight;
    const scalar cosAlpha = sign((npAlpha - sp) & dv) * mag(npAlpha - sp);
    const scalar alpha =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosAlpha / (currThickness + VSMALL))
            )
        );

    //- calculate the height of the layer at the neighbour
    //- point
    const point npBeta = help::nearestPointOnTheEdge(ep, sp, points[nhe[1]]);
    const scalar currNeiHeight = mag(npBeta - points[nhe[1]]);
    scalar suggestedNeiHeight = currNeiHeight;
    const scalar cosBeta = sign((npBeta - ep) & -dv) * mag(npBeta - ep);
    const scalar beta =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosBeta / (currNeiThickness + VSMALL))
            )
        );

    //- check if the current thickness is Ok for the local curvature
    if( (alpha + beta) < M_PI )
    {
        const scalar gamma = M_PI - (alpha + beta);
        const scalar sinGamma = Foam::max(SMALL, Foam::sin(gamma));
        const scalar sinAlpha = Foam::max(SMALL, Foam::sin(alpha));
        const scalar sinBeta = Foam::max(SMALL, Foam::sin(beta));

        //- max allowed thickness and layer height due to curvature
        retThickness =
            Foam::min
            (
                retThickness,
                featureSizeFactor_ * magDv * sinBeta / sinGamma
            );
        retHeight *= (retThickness / (currThickness + VSMALL));

        //- max allowed neighbour hair thickness
        //- and layer height due to curvature
        suggestedNeiThickness =
            Foam::min
            (
                suggestedNeiThickness,
                featureSizeFactor_ * magDv * sinAlpha / sinGamma
            );
        suggestedNeiHeight *=
            (suggestedNeiThickness / (currNeiThickness + VSMALL));
    }

    //- check the height variation
    const scalar tanVal = (retHeight - suggestedNeiHeight) / (magDv + VSMALL);

    if( tanVal > relThicknessTol_ )
    {
        retHeight = suggestedNeiHeight + relThicknessTol_ * magDv;

        retThickness = (retHeight / currHeight) * currThickness;
    }

    return retThickness;
}

scalar boundaryLayerOptimisation::calculateThicknessOverCell
(
    const label heI,
    const label cellI,
    const label baseFaceI
) const
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

    const face& bf = faces[baseFaceI];

    const edge& he = hairEdges_[heI];

    const point& sp = points[he[0]];
    const point& ep = points[he[1]];

    scalar maxThickness = he.mag(points);

    //- the base face must not contain the hair edge
    //- this is the case at exitting layers
    forAll(bf, eI)
        if( bf.faceEdge(eI) == he )
            return maxThickness;

    forAll(c, fI)
    {
        if( c[fI] == baseFaceI )
            continue;

        const face& f = faces[c[fI]];

        if( help::shareAnEdge(bf, f) && (f.which(he.start()) == -1) )
        {
            point intersection;

            if( !help::lineFaceIntersection(sp, ep, f, points, intersection) )
                continue;

            const scalar maxDist = featureSizeFactor_ * mag(intersection - sp);

            maxThickness =
                Foam::min(maxThickness, maxDist);
        }
    }

    return maxThickness;
}

void boundaryLayerOptimisation::optimiseThicknessVariation
(
    const direction edgeType
)
{
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const label start = mesh_.nInternalFaces();
    const labelList& faceOwner = mse.faceOwners();

    vectorField hairDirections(hairEdges_.size());
    scalarField hairLength(hairEdges_.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        vector n = hairEdges_[hairEdgeI].vec(points);

        hairLength[hairEdgeI] = (Foam::mag(n) + VSMALL);
        hairDirections[hairEdgeI] = n / hairLength[hairEdgeI];
    }

    //- reduce thickness of the layer
    //- such that the variation of layer thickness
    //- It is an iterative process where the layer is thinned in the regions
    //- where the tangent is greater than the tolerance value or the curvature
    //- permits thicker boundary layers.
    boolList activeHairEdge(hairEdges_.size(), true);
    bool changed;
    label nIter(0);
    do
    {
        changed = false;

        boolList modifiedHairEdge(hairEdges_.size(), false);
        boolList influencedEdges(hairEdges_.size(), false);

        //- check if the hair edge intersects some other face in the cells
        //- attached to the hair edge
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            if
            (
                (hairEdgeType_[hairEdgeI] & edgeType) &&
                activeHairEdge[hairEdgeI]
            )
            {
                const edge& he = hairEdges_[hairEdgeI];

                const label bpI = bp[he.start()];

                scalar maxThickness = hairLength[hairEdgeI];

                DynList<label, 64> influencers;

                forAllRow(pFaces, bpI, pfI)
                {
                    const label bfI = pFaces(bpI, pfI);

                    const face& bf = bFaces[bfI];
                    if( bf.which(he.end()) >= 0 )
                        continue;

                    const label baseFaceI = start + bfI;
                    const label cOwn = faceOwner[bfI];

                    //- check if there exist any self-intersections
                    maxThickness =
                        Foam::min
                        (
                            maxThickness,
                            calculateThicknessOverCell
                            (
                                hairEdgeI,
                                cOwn,
                                baseFaceI
                            )
                        );

                    //- check thickness variation over all hair edges
                    DynList<edge> hairEdges;
                    hairEdgesAtBndFace(faceOwner[bfI], baseFaceI, hairEdges);

                    forAll(bf, pI)
                    {
                        const edge& nhe = hairEdges[pI];

                        const label bpJ = bp[nhe.start()];

                        forAllRow(hairEdgesAtBndPoint_, bpJ, peJ)
                        {
                            const label heJ = hairEdgesAtBndPoint_(bpJ, peJ);

                            if( hairEdgeI == heJ )
                                continue;

                            if( nhe == hairEdges_[heJ] )
                            {
                                influencers.append(heJ);

                                const scalar edgeThickness =
                                    calculateThickness(hairEdgeI, heJ);

                                maxThickness =
                                    Foam::min(maxThickness, edgeThickness);
                            }
                        }
                    }
                }

                if( hairLength[hairEdgeI] > maxThickness )
                {
                    //- make the hair edge shorter
                    hairLength[hairEdgeI] = maxThickness;
                    modifiedHairEdge[hairEdgeI] = true;
                    changed = true;

                    thinnedHairEdge_[hairEdgeI] = true;

                    forAll(influencers, i)
                        influencedEdges[influencers[i]] = true;
                }
            }
        }
/*
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            const scalar magN = hairLength[hairEdgeI];

            if( magN < VSMALL )
                FatalErrorIn
                (
                    "void boundaryLayerOptimisation::optimiseThicknessVariation"
                    "(const direction, const label, const scalar, const scalar)"
                ) << "Zero layer thickness at hair edge " << hairEdgeI
                  << ". Exitting..." << exit(FatalError);

            if( hairEdgeType_[hairEdgeI] & edgeType )
            {
                forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                {
                    const label hairEdgeJ =
                        hairEdgesNearHairEdge_(hairEdgeI, nheI);

                    if( !activeHairEdge[hairEdgeJ] )
                        continue;

                    const scalar maxThickness =
                        calculateThickness
                        (
                            hairEdgeI,
                            hairEdgeJ
                        );

                    if( hairLength[hairEdgeI] > maxThickness )
                    {
                        //- make the hair edge shorter
                        hairLength[hairEdgeI] = maxThickness;

                        changed = true;
                        thinnedHairEdge_[hairEdgeI] = true;
                        modifiedHairEdge[hairEdgeI] = true;
                    }
                }
            }
        }
        */

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();

            const edgeList& edges = mesh_.addressingData().edges();
            const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
            const VRWGraph& edgesAtProcs =
                mesh_.addressingData().edgeAtProcs();
            const labelLongList& globalEdgeLabel =
                mesh_.addressingData().globalEdgeLabel();
            const Map<label>& globalToLocalEdge =
                mesh_.addressingData().globalToLocalEdgeAddressing();
            const DynList<label>& eNeiProcs =
                mesh_.addressingData().edgeNeiProcs();

            std::map<label, LongList<labelledScalar> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                    if( !(hairEdgeType_[hairEdgeI] & edgeType) )
                        continue;

                    const edge& he = hairEdges_[hairEdgeI];

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            labelledScalar lScalar
                            (
                                globalEdgeLabel[edgeI],
                                hairLength[hairEdgeI]
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<labelledScalar>& dts =
                                    exchangeData[neiProc];

                                dts.append(lScalar);
                            }
                        }
                    }
                }
            }

            LongList<labelledScalar> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelledScalar& lScalar = receivedData[i];
                const label edgeI = globalToLocalEdge[lScalar.scalarLabel()];
                const edge& e = edges[edgeI];

                bool found(false);
                for(label pI=0;pI<2;++pI)
                {
                    const label bpI = bp[e[pI]];

                    if( bpI < 0 )
                        continue;

                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( hairEdges_[hairEdgeI] == e )
                        {
                            if( lScalar.value() < hairLength[hairEdgeI] )
                            {
                                hairLength[hairEdgeI] = lScalar.value();
                                changed = true;
                                thinnedHairEdge_[hairEdgeI] = true;
                                modifiedHairEdge[hairEdgeI] = true;
                            }

                            found = true;
                            break;
                        }
                    }

                    if( found )
                        break;
                }
            }
        }

        //- reduce the information over all processors
        reduce(changed, maxOp<bool>());

        if( !changed )
            break;

        //- move boundary vertices to the new positions
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            if( hairEdgeType_[hairEdgeI] & edgeType )
            {
                const edge& he = hairEdges_[hairEdgeI];
                const vector& hv = hairDirections[hairEdgeI];
                const point& s = points[he.start()];

                points[he.end()] = s + hairLength[hairEdgeI] * hv;
            }
        }

        //- mark edges which may be changed
        activeHairEdge.transfer(influencedEdges);
    } while( changed && (++nIter < 1000) );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
