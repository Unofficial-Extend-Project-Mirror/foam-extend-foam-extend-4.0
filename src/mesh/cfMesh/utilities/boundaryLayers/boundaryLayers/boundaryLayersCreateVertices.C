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

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "labelledScalar.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::findPatchVertices
(
    const boolList& treatPatches,
    List<direction>& pVertices
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    pVertices.setSize(pPatches.size());
    pVertices = NONE;

    # ifdef USE_OMP
    # pragma omp parallel for if( pPatches.size() > 1000 ) \
    schedule(dynamic, Foam::max(10, pPatches.size()/(2*omp_get_num_threads())))
    # endif
    forAll(pPatches, bpI)
    {
        bool hasTreated(false);
        bool hasNotTreated(false);

        forAllRow(pPatches, bpI, patchI)
        {
            const label patch = pPatches(bpI, patchI);
            if( treatPatches[patch] )
            {
                hasTreated = true;
            }
            else
            {
                hasNotTreated = true;
            }
        }

        if( hasTreated )
        {
            pVertices[bpI] |= PATCHNODE;

            if( hasNotTreated )
                pVertices[bpI] |= EDGENODE;
        }
    }

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        forAll(pVertices, bpI)
            if( pVertices[bpI] && (bpAtProcs.sizeOfRow(bpI) != 0) )
                pVertices[bpI] |= PARALLELBOUNDARY;
    }
}

point boundaryLayers::createNewVertex
(
    const label bpI,
    const boolList& treatPatches,
    const List<direction>& patchVertex
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const VRWGraph& pointPoints = mse.pointPoints();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = mesh_.points();

    # ifdef DEBUGLayer
    Info << "Creating new vertex for boundary vertex " << bpI << endl;
    Info << "Global vertex label " << bPoints[bpI] << endl;
    # endif

    vector normal(vector::zero);
    scalar dist(VGREAT);
    const point& p = points[bPoints[bpI]];
    if( patchVertex[bpI] & EDGENODE )
    {
        # ifdef DEBUGLayer
        Info << "Vertex is on the border" << endl;
        # endif

        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn
                (
                    pPatches(bpI, patchI)
                );

        if( otherPatches.size() == 1 )
        {
            //- vertex is on an edge
            # ifdef DEBUGLayer
            Info << "Vertex is on an edge" << endl;
            # endif
            vector v(vector::zero);

            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if( treatPatches[patchLabel] )
                {
                    normal += f.normal(points);
                }
                else
                {
                    v += f.normal(points);
                }
            }

            const scalar magV = mag(v) + VSMALL;
            v /= magV;

            normal -= (normal & v) * v;

            const scalar magN = mag(normal) + VSMALL;
            normal /= magN;

            forAllRow(pointPoints, bpI, ppI)
            {
                if( patchVertex[pointPoints(bpI, ppI)] )
                    continue;

                const vector vec = points[bPoints[pointPoints(bpI, ppI)]] - p;
                const scalar prod = 0.5 * mag(vec & normal);

                if( prod < dist )
                    dist = prod;
            }
        }
        else if( otherPatches.size() == 2 )
        {
            # ifdef DEBUGLayer
            Info << "Vertex is a corner" << endl;
            # endif

            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                    if( !pPatches.contains(bpJ, otherPatches[opI]) )
                    {
                        found = false;
                        break;
                    }

                if( found )
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if( otherVertex == -1 )
            {
                FatalErrorIn
                (
                    "void boundaryLayers::createNewVertices"
                    "("
                        "const boolList& treatPatches,"
                        "labelList& newLabelForVertex"
                    ")"
                ) << "Cannot find moving vertex!" << exit(FatalError);
            }

            //- normal vector is co-linear with that edge
            normal = p - points[bPoints[otherVertex]];
            dist = 0.5 * mag(normal) + VSMALL;

            normal /= 2.0 * dist;
        }
        else
        {
            FatalErrorIn
            (
                "void boundaryLayers::createNewVertices"
                "("
                    "const boolList& treatPatches,"
                    "labelList& newLabelForVertex"
                ") const"
            ) << "There are more than 3 patches meeting at this vertex!"
                << pPatches[bpI] << abort(FatalError);
        }

        //- limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if( otherPatches.contains(boundaryFacePatches[faceLabel]) )
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if( pos != -1 )
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if( dst < dist )
                        dist = 0.9 * dst;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayers::createNewVertices"
                        "("
                            "const boolList& treatPatches,"
                            "labelList& newLabelForVertex"
                        ") const"
                    ) << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        normal = pNormals[bpI];

        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
            0.5 * mag
            (
                points[bPoints[pointPoints(bpI, ppI)]] -
                p
            );

            if( d < dist )
                dist = d;
        }
    }

    //- create new vertex
    # ifdef DEBUGLayer
    Info << "Normal for vertex " << bpI << " is " << normal << endl;
    Info << "Distance is " << dist << endl;
    # endif

    dist = Foam::max(dist, VSMALL);

    const point newP = p - dist * normal;

    if( help::isnan(newP) || help::isinf(newP) )
        return p;

    return newP;
}

void boundaryLayers::createNewVertices(const boolList& treatPatches)
{
    Info << "Creating vertices for layer cells" << endl;

    List<direction> patchVertex;
    findPatchVertices(treatPatches, patchVertex);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    //- the following is needed for parallel runs
    //- it is ugly, but must stay for now :(
    if( Pstream::parRun() )
    {
        mse.pointNormals();
        mse.pointPoints();
    }

    pointFieldPMG& points = mesh_.points();

    label nExtrudedVertices(0);
    forAll(patchVertex, bpI)
        if( patchVertex[bpI] )
            ++nExtrudedVertices;

    points.setSize(points.size() + nExtrudedVertices);

    labelLongList procPoints;
    forAll(bPoints, bpI)
        if( patchVertex[bpI] )
        {
            if( patchVertex[bpI] & PARALLELBOUNDARY )
            {
                procPoints.append(bpI);
                continue;
            }

            points[nPoints_] = createNewVertex(bpI, treatPatches, patchVertex);
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
            ++nPoints_;
        }

    if( Pstream::parRun() )
    {
        createNewPartitionVerticesParallel
        (
            procPoints,
            patchVertex,
            treatPatches
        );

        createNewEdgeVerticesParallel
        (
            procPoints,
            patchVertex,
            treatPatches
        );
    }

    //- swap coordinates of new and old points
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];
        if( pLabel != -1 )
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }

    if( nPoints_ != points.size() )
        FatalErrorIn
        (
            "void boundaryLayers::createNewVertices("
            "const meshSurfaceEngine& mse,"
            "const boolList& treatPatches,"
            "labelList& newLabelForVertex)"
        ) << "Number of vertices " << nPoints_
            << " does not match the list size "
            << abort(FatalError);

    Info << "Finished creating layer vertices" << endl;
}

void boundaryLayers::createNewVertices(const labelList& patchLabels)
{
    otherVrts_.clear();

    patchKey_.setSize(mesh_.boundaries().size());
    patchKey_ = -1;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    //- the following is needed for parallel runs
    //- it is ugly, but must stay for now :(
    mse.boundaryFaces();
    mse.pointNormals();
    mse.pointFaces();
    mse.pointPoints();

    pointFieldPMG& points = mesh_.points();
    boolList treatPatches(mesh_.boundaries().size());
    List<direction> patchVertex(bPoints.size());

    //- make sure than the points are never re-allocated during the process
    points.reserve(points.size() + 2 * bPoints.size());

    //- generate new layer vertices for each patch
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        treatPatches = false;

        bool treat(true);
        forAll(treatPatchesWithPatch_[pLabel], pI)
        {
            const label otherPatch = treatPatchesWithPatch_[pLabel][pI];
            treatPatches[otherPatch] = true;

            if( patchKey_[otherPatch] == -1 )
            {
                patchKey_[otherPatch] = patchI;
            }
            else
            {
                treat = false;
            }
        }

        if( !treat )
            continue;

        const label pKey = patchKey_[pLabel];

        //- classify vertices belonging to this patch
        findPatchVertices(treatPatches, patchVertex);

        //- create indices and allocate maps for new points
        labelLongList procPoints, patchPoints;
        forAll(bPoints, bpI)
        {
            if( !patchVertex[bpI] )
                continue;

            //- skip vertices at parallel boundaries
            if( patchVertex[bpI] & PARALLELBOUNDARY )
            {
                procPoints.append(bpI);

                continue;
            }

            patchPoints.append(bpI);
            const label pointI = bPoints[bpI];

            if( patchVertex[bpI] & EDGENODE )
            {
                if( otherVrts_.find(pointI) == otherVrts_.end() )
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                otherVrts_[pointI].insert(std::make_pair(pr, nPoints_++));
            }
            else
            {
                //- this the only new point
                newLabelForVertex_[pointI] = nPoints_++;
            }
        }

        //- set the size of points
        points.setSize(nPoints_);

        //- calculate coordinates of new points
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(patchPoints, i)
        {
            const label bpI = patchPoints[i];

            const label pointI = bPoints[bpI];

            //- create new point
            const point p = createNewVertex(bpI, treatPatches, patchVertex);

            if( patchVertex[bpI] & EDGENODE )
            {
                //- set the new point or an edge point
                if( otherVrts_.find(pointI) == otherVrts_.end() )
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                const label npI = otherVrts_[pointI][pr];
                points[npI] = p;
            }
            else
            {
                //- set the new point
                points[newLabelForVertex_[pointI]] = p;
            }
        }

        if( Pstream::parRun() )
        {
            points.setSize(nPoints_+procPoints.size());

            createNewPartitionVerticesParallel
            (
                procPoints,
                patchVertex,
                treatPatches
            );

            createNewEdgeVerticesParallel
            (
                procPoints,
                patchVertex,
                treatPatches
            );
        }
    }

    //- create missing vertices for edge and corner vertices
    //- they should be stored in the otherNodes map
    forAll(bPoints, bpI)
    {
        const label pointI = bPoints[bpI];

        if( otherVrts_.find(pointI) == otherVrts_.end() )
            continue;

        const point& p = points[pointI];

        std::map<std::pair<label, label>, label>& m = otherVrts_[pointI];
        DynList<label> usedPatches;
        DynList<label> newNodeLabel;
        DynList<vector> newPatchPenetrationVector;
        forAllRow(pPatches, bpI, patchI)
        {
            const label pKey = patchKey_[pPatches(bpI, patchI)];
            const std::pair<label, label> pr(pKey, pKey);
            const std::map<std::pair<label, label>, label>::const_iterator it =
                m.find(pr);
            if( (it != m.end()) && !usedPatches.contains(pKey) )
            {
                usedPatches.append(pKey);
                newNodeLabel.append(it->second);
                newPatchPenetrationVector.append(points[it->second] - p);
            }
        }

        if( newNodeLabel.size() == 1 )
        {
            //- only one patch is treated
            newLabelForVertex_[pointI] = newNodeLabel[0];
            otherVrts_.erase(pointI);
        }
        else if( newNodeLabel.size() == 2 )
        {
            //- point is located at an extrusion edge
            //- create the new position for the existing point
            point newP(p);
            newP += newPatchPenetrationVector[0];
            newP += newPatchPenetrationVector[1];

            if( !help::isnan(newP) && !help::isinf(newP) )
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }
            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else if( newNodeLabel.size() == 3 )
        {
            //- point is located at an extrusion corner
            //- create 3 points and the new position for the existing point
            point newP(p);
            for(label i=0;i<3;++i)
            {
                newP += newPatchPenetrationVector[i];
                for(label j=i+1;j<3;++j)
                {
                    const point np =
                        p + newPatchPenetrationVector[i] +
                        newPatchPenetrationVector[j];

                    if( !help::isnan(np) && !help::isinf(np) )
                    {
                        points.append(np);
                    }
                    else
                    {
                        points.append(p);
                    }

                    m.insert
                    (
                        std::make_pair
                        (
                            std::make_pair(usedPatches[i], usedPatches[j]),
                            nPoints_
                        )
                    );
                    ++nPoints_;
                }
            }

            //- create new position for the existing point
            if( !help::isnan(newP) && !help::isinf(newP) )
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }

            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else
        {
            FatalErrorIn
            (
                "void boundaryLayers::createNewVertices("
                "const labelList& patchLabels, labelLongList& newLabelForVertex,"
                "std::map<label, std::map<std::pair<label, label>, label> >&)"
            ) << "Boundary node " << bpI << " is not at an edge!"
                << abort(FatalError);
        }
    }

    //- swap coordinates of new and old points
    # ifdef USE_OMP
    # pragma omp parallel for if( bPoints.size() > 1000 ) \
    schedule(dynamic, 100)
    # endif
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];

        if( pLabel != -1 )
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }
}

void boundaryLayers::createNewPartitionVerticesParallel
(
    const labelLongList& procPoints,
    const List<direction>& pVertices,
    const boolList& /*treatPatches*/
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    scalarField penetrationDistances(bPoints.size(), VGREAT);

    std::map<label, LongList<labelledScalar> > exchangeDistances;

    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeDistances.find(neiProc) == exchangeDistances.end() )
            {
                exchangeDistances.insert
                (
                    std::make_pair(neiProc, LongList<labelledScalar>())
                );
            }
        }

        if( pVertices[bpI] & EDGENODE )
            continue;

        scalar dist(VGREAT);
        const point& p = points[bPoints[bpI]];
        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
                0.5 * mag(points[bPoints[pointPoints(bpI, ppI)]] - p);

            if( d < dist )
                dist = d;
        }

        penetrationDistances[bpI] = dist;

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeDistances[neiProc].append
            (
                labelledScalar(globalPointLabel[bpI], dist)
            );
        }
    }

    //- exchange distances with other processors
    LongList<labelledScalar> receivedData;
    help::exchangeMap(exchangeDistances, receivedData);
    forAll(receivedData, i)
    {
        const label bpI = globalToLocal[receivedData[i].scalarLabel()];

        if( penetrationDistances[bpI] > receivedData[i].value() )
            penetrationDistances[bpI] = receivedData[i].value();
    }

    //- Finally, create the points
    const vectorField& pNormals = mse.pointNormals();
    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];

        if( pVertices[bpI] & EDGENODE )
            continue;

        const point& p = points[bPoints[bpI]];
        const point np = p - pNormals[bpI] * penetrationDistances[bpI];
        if( !help::isnan(np) && !help::isinf(np) )
        {
            points[nPoints_] = np;
        }
        else
        {
            points[nPoints_] = p;
        }
        newLabelForVertex_[bPoints[bpI]] = nPoints_;
        ++nPoints_;
    }
}

void boundaryLayers::createNewEdgeVerticesParallel
(
    const labelLongList& procPoints,
    const List<direction>& pVertices,
    const boolList& treatPatches
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    DynList<label> neiProcs;
    labelLongList edgePoints;
    Map<label> bpToEdgePoint;
    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            neiProcs.appendIfNotIn(neiProc);
        }

        if( pVertices[bpI] & EDGENODE )
        {
            bpToEdgePoint.insert(bpI, edgePoints.size());
            edgePoints.append(bpI);
        }
    }

    if( returnReduce(edgePoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    scalarField dist(edgePoints.size(), VGREAT);
    vectorField normal(edgePoints.size(), vector::zero);
    vectorField v(edgePoints.size(), vector::zero);

    label pKey(-1);
    if( patchKey_.size() )
    {
        forAll(treatPatches, patchI)
            if( treatPatches[patchI] )
            {
                pKey = patchKey_[patchI];
                break;
            }
    }

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        const point& p = points[bPoints[bpI]];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        //- find local values of normals and v
        if( otherPatches.size() == 1 )
        {
            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if( treatPatches[patchLabel] )
                {
                    normal[epI] += f.normal(points);
                }
                else
                {
                    v[epI] += f.normal(points);
                }
            }
        }
        else if( otherPatches.size() == 2 )
        {
            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                    if( !pPatches.contains(bpJ, otherPatches[opI]) )
                    {
                        found = false;
                        break;
                    }

                if( found )
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if( otherVertex == -1 )
                continue;

            //- normal vector is co-linear with that edge
            normal[epI] = p - points[bPoints[otherVertex]];
            dist[epI] = mag(normal[epI]);
        }
        else
        {
            FatalErrorIn
            (
                "void boundaryLayers::createNewEdgeVerticesParallel("
                    "const labelLongList& procPoints,"
                    "const List<direction>& pVertices,"
                    "const boolList& treatPatches,"
                    "labelList& newLabelForVertex"
                ") const"
            ) << "There are more than 3 patches meeting at this vertex!"
                << abort(FatalError);
        }
    }

    //- prepare normals and v for sending to other procs
    std::map<label, LongList<labelledPoint> > exchangeNormals;
    forAll(neiProcs, procI)
        exchangeNormals.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledPoint>())
        );

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- store values in the list for sending
            LongList<labelledPoint>& dataToSend = exchangeNormals[neiProc];
            dataToSend.append
            (
                labelledPoint(globalPointLabel[bpI], normal[epI])
            );
            dataToSend.append(labelledPoint(globalPointLabel[bpI], v[epI]));
        }
    }

    //- exchange data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeNormals, receivedData);
    exchangeNormals.clear();

    label counter(0);
    while( counter < receivedData.size() )
    {
        const labelledPoint& otherNormal = receivedData[counter++];
        const labelledPoint& otherV = receivedData[counter++];

        const label bpI = globalToLocal[otherNormal.pointLabel()];
        normal[bpToEdgePoint[bpI]] += otherNormal.coordinates();
        v[bpToEdgePoint[bpI]] += otherV.coordinates();
    }

    //- calculate normals
    forAll(normal, epI)
    {
        const label bpI = edgePoints[epI];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        if( otherPatches.size() == 1 )
        {
            const scalar magV = mag(v[epI]) + VSMALL;
            v[epI] /= magV;
            normal[epI] -= (normal[epI] & v[epI]) * v[epI];
        }

        const scalar magN = mag(normal[epI]) + VSMALL;
        normal[epI] /= magN;
    }

    //- calculate distances
    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        const point& p = points[bPoints[bpI]];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        if( otherPatches.size() == 1 )
        {
            forAllRow(pointPoints, bpI, ppI)
            {
                if( pVertices[pointPoints(bpI, ppI)] )
                    continue;

                const vector vec = points[bPoints[pointPoints(bpI, ppI)]] - p;
                const scalar prod = 0.5 * mag(vec & normal[epI]);

                if( prod < dist[epI] )
                    dist[epI] = prod;
            }
        }

        //- limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if( otherPatches.contains(boundaryFacePatches[faceLabel]) )
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if( pos != -1 )
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if( dst < dist[epI] )
                        dist[epI] = 0.9 * dst;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayers::createNewEdgeVerticesParallel"
                        "("
                            "const labelLongList& procPoints,"
                            "const List<direction>& pVertices,"
                            "const boolList& treatPatches,"
                            "labelList& newLabelForVertex"
                        ") const"
                    ) << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }

    //- exchange distances with other processors
    std::map<label, LongList<labelledScalar> > exchangeDistances;
    forAll(neiProcs, procI)
        exchangeDistances.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledScalar>())
        );

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            LongList<labelledScalar>& ls = exchangeDistances[neiProc];
            ls.append(labelledScalar(globalPointLabel[bpI], dist[epI]));
        }
    }

    //- exchange distances with other processors
    LongList<labelledScalar> receivedDistances;
    help::exchangeMap(exchangeDistances, receivedDistances);
    exchangeDistances.clear();

    forAll(receivedDistances, i)
    {
        const label bpI = globalToLocal[receivedDistances[i].scalarLabel()];
        const label epI = bpToEdgePoint[bpI];
        if( dist[epI] > receivedDistances[i].value() )
            dist[epI] = receivedDistances[i].value();
    }

    //- Finally, create new points
    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];

        const point& p = points[bPoints[bpI]];
        const point np = p - normal[epI] * dist[epI];
        if( !help::isnan(np) && !help::isinf(np) )
        {
            points[nPoints_] = np;
        }
        else
        {
            points[nPoints_] = p;
        }

        if( pKey == -1 )
        {
            //- extrusion for one patch in a single go
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
        }
        else
        {
            const label pointI = bPoints[bpI];

            if( otherVrts_.find(pointI) == otherVrts_.end() )
            {
                std::map<std::pair<label, label>, label> m;
                otherVrts_.insert(std::make_pair(pointI, m));
            }

            std::pair<label, label> pr(pKey, pKey);
            otherVrts_[pointI].insert(std::make_pair(pr, nPoints_));
        }
        ++nPoints_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
