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

#include "meshSurfaceCheckEdgeTypes.H"
#include "meshSurfaceEngine.H"
#include "boolList.H"
#include "demandDrivenData.H"
#include "helperFunctionsPar.H"
#include "triangle.H"
#include "tetrahedron.H"
#include "labelledPoint.H"

#include <map>
# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCheckEdgeTypes::classifyEdges()
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bp = surfaceEngine_.bp();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const edgeList& edges = surfaceEngine_.edges();
    const VRWGraph& edgeFaces = surfaceEngine_.edgeFaces();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const labelList& facePatch = surfaceEngine_.boundaryFacePatches();
    const vectorField& pNormals = surfaceEngine_.pointNormals();
    const vectorField& fCentres = surfaceEngine_.faceCentres();
    const vectorField& fNormals = surfaceEngine_.faceNormals();

    boolList problematicPoint(pointFaces.size());
    edgeTypes_.setSize(edges.size());

    # ifdef USE_OMP
    label nThreads = 3 * omp_get_num_procs();
    if( bFaces.size() < 1000 )
        nThreads = 1;
    # endif

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(edgeTypes_, edgeI)
            edgeTypes_[edgeI] = NONE;

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(problematicPoint, pointI)
            problematicPoint[pointI] = false;

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp for schedule(static, 1)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            //- check if point normal are consistent with face normals
            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];
                const label bpNext = bp[bf.nextLabel(pI)];

                if( (pNormals[bpI] & fNormals[bfI]) < VSMALL )
                    problematicPoint[bpI] = true;

                //- contruct a triangle from a face edge and centre point
                const triangle<point, point> tria
                (
                    points[bf[pI]],
                    points[bf.nextLabel(pI)],
                    fCentres[bfI]
                );

                //- check if the normal of the triangle is consistent with
                //- the face normal
                if( (tria.normal() & fNormals[bfI]) < VSMALL )
                {
                    problematicPoint[bpI] = true;
                    problematicPoint[bpNext] = true;
                }

                //- check if the face normal is consistent with point normals
                if( (tria.normal() & pNormals[bpI]) < VSMALL )
                    problematicPoint[bpI] = true;

                if( (tria.normal() & pNormals[bpNext]) < VSMALL )
                    problematicPoint[bpNext] = true;
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier

        //- start checking feature edges
        # pragma omp for schedule(static, 1)
        # endif
        forAll(edgeFaces, edgeI)
        {
            if( edgeFaces.sizeOfRow(edgeI) == 2 )
            {
                const label f0 = edgeFaces(edgeI, 0);
                const label f1 = edgeFaces(edgeI, 1);

                if( facePatch[f0] == facePatch[f1] )
                    edgeTypes_[edgeI] |= PATCHEDGE;

                const edge e = edges[edgeI];

                //- check if the surface is tangled there
                if( problematicPoint[bp[e.start()]] )
                {
                    edgeTypes_[edgeI] |= UNDETERMINED;
                    continue;
                }

                if( problematicPoint[bp[e.end()]] )
                {
                    edgeTypes_[edgeI] |= UNDETERMINED;
                    continue;
                }

                //- check the volumes pof tets which can be formed at the edge
                const tetrahedron<point, point> tet0
                (
                    points[e.start()],
                    points[e.end()],
                    fCentres[f0],
                    fCentres[f1]
                );

                if( tet0.mag() > -VSMALL )
                {
                    edgeTypes_[edgeI] |= CONCAVEEDGE;
                    continue;
                }

                const tetrahedron<point, point> tet1
                (
                    points[e.end()],
                    points[e.start()],
                    fCentres[f1],
                    fCentres[f0]
                );

                if( tet1.mag() > -VSMALL )
                {
                    edgeTypes_[edgeI] |= CONCAVEEDGE;
                    continue;
                }

                edgeTypes_[edgeI] |= CONVEXEDGE;
            }
        }
    }

    if( Pstream::parRun() )
    {
        const labelList& globalPointLabel =
            surfaceEngine_.globalBoundaryPointLabel();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const DynList<label>& bpNeiProcs = surfaceEngine_.bpNeiProcs();

        //- make sure that problematic points
        //- are consistent ove processor boundaries
        std::map<label, labelLongList> exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData.insert
            (
                std::make_pair(bpNeiProcs[i], labelLongList())
            );

        forAllConstIter(Map<label>, globalToLocal, bpIter)
        {
            const label bpI = bpIter();

            if( !problematicPoint[bpI] )
                continue;

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProcs = bpAtProcs(bpI, i);

                if( neiProcs == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProcs].append(globalPointLabel[bpI]);
            }
        }

        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        forAll(receiveData, i)
            problematicPoint[globalToLocal[receiveData[i]]] = true;

        //- check if the edge at processor boudandaries are concave or convex
        const labelList& globalEdgeLabel =
            surfaceEngine_.globalBoundaryEdgeLabel();
        const Map<label>& otherProc = surfaceEngine_.otherEdgeFaceAtProc();
        const Map<label>& otherPatch = surfaceEngine_.otherEdgeFacePatch();
        const Map<label>& globalToLocalEdge =
            surfaceEngine_.globalToLocalBndEdgeAddressing();

        std::map<label, LongList<labelledPoint> > exchangeFaceCentres;
        forAll(surfaceEngine_.beNeiProcs(), i)
        {
            const label neiProc = surfaceEngine_.beNeiProcs()[i];

            exchangeFaceCentres.insert
            (
                std::make_pair(neiProc, LongList<labelledPoint>())
            );
        }

        forAllConstIter(Map<label>, otherPatch, eIter)
        {
            if( eIter() == facePatch[edgeFaces(eIter.key(), 0)] )
                continue;

            const edge& e = edges[eIter.key()];
            if( problematicPoint[e.start()] || problematicPoint[e.end()] )
            {
                edgeTypes_[eIter.key()] |= UNDETERMINED;
                continue;
            }

            const label neiProcs = otherProc[eIter.key()];
            exchangeFaceCentres[neiProcs].append
            (
                labelledPoint
                (
                    globalEdgeLabel[eIter.key()],
                    fCentres[edgeFaces(eIter.key(), 0)]
                )
            );
        }

        LongList<labelledPoint> receiveCentres;
        help::exchangeMap(exchangeFaceCentres, receiveCentres);

        forAll(receiveCentres, i)
        {
            const labelledPoint& lp = receiveCentres[i];
            const label edgeI = globalToLocalEdge[lp.pointLabel()];
            const vector fCentre = lp.coordinates();

            const edge& e = edges[edgeI];
            const label f0 = edgeFaces(edgeI, 0);

            //- check the volumes pof tets
            //- which can be formed at the edge
            tetrahedron<point, point> tet0
            (
                points[e.start()],
                points[e.end()],
                fCentres[f0],
                fCentre
            );

            if ( tet0.mag() > -VSMALL )
            {
                edgeTypes_[edgeI] |= CONCAVEEDGE;
                continue;
            }

            tetrahedron<point, point> tet1
            (
                points[e.end()],
                points[e.start()],
                fCentre,
                fCentres[f0]
            );

            if ( tet1.mag() > -VSMALL )
            {
                edgeTypes_[edgeI] |= CONCAVEEDGE;
                continue;
            }

            edgeTypes_[edgeI] |= CONVEXEDGE;
        }
    }

    # ifdef DEBUGClassifyEdges
    polyMeshGen& mesh_ = const_cast<polyMeshGen&>(surfaceEngine_.mesh());
    const label badVertices = mesh_.addPointSubset("invertedVertices");
    forAll(problematicPoint, bpI)
        if( problematicPoint[bpI] )
            mesh_.addPointToSubset
            (
                badVertices,
                surfaceEngine_.boundaryPoints()[bpI]
            );

    forAll(edgeTypes_, edgeI)
    {
        if( edgeTypes_[edgeI] & CONVEXEDGE )
        {
            Info <<"Edge " << edgeI << " is convex" << endl;
        }
        else if( edgeTypes_[edgeI] & CONCAVEEDGE )
        {
            Info << "Edge " << edgeI << " is concave" << endl;
        }
        else if( edgeTypes_[edgeI] & UNDETERMINED )
        {
            Info << "Edge " << edgeI << " is not determined" << endl;
        }
        else if( edgeTypes_[edgeI] & PATCHEDGE )
        {
            Info << "Edge " << edgeI << " is a patch edge" << endl;
        }
        else
        {
            Info << "Drekec spekec" << edgeI << endl;
        }
    }
    # endif
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceCheckEdgeTypes::meshSurfaceCheckEdgeTypes
(
    const meshSurfaceEngine& mse
)
:
    surfaceEngine_(mse),
    edgeTypes_()
{
    classifyEdges();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceCheckEdgeTypes::~meshSurfaceCheckEdgeTypes()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCheckEdgeTypes::convexEdges(labelLongList& convexEdges) const
{
    convexEdges.clear();

    forAll(edgeTypes_, eI)
    {
        if( edgeTypes_[eI] & CONVEXEDGE )
            convexEdges.append(eI);
    }
}

void meshSurfaceCheckEdgeTypes::concaveEdges(labelLongList& concaveEdges) const
{
    concaveEdges.clear();

    forAll(edgeTypes_, eI)
    {
        if( edgeTypes_[eI] & CONCAVEEDGE )
            concaveEdges.append(eI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
