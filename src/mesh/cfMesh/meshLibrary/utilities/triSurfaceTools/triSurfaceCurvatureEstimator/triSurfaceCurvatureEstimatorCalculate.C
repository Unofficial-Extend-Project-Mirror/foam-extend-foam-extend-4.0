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

#include "triSurfaceCurvatureEstimator.H"
#include "matrix3D.H"
#include "quadricFitting.H"
#include "HashSet.H"
#include "boolList.H"
#include "Map.H"
#include "DynList.H"

#include "OFstream.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCurvatureEstimator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeSurfaceToVTK
(
    OFstream& file,
    const triSurf& surf
)
{
    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    const pointField& points = surf.points();
    file << "POINTS " << points.size() << " float\n";
    forAll(points, pI)
    {
        const point& p = points[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';

        if( pI % 5 == 0 )
            file << "\n";
    }

    //- write triangles
    file << "\n";
    file << "\nPOLYGONS " << surf.size() << " " << 4*surf.size() << endl;
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];
        file << 3 << " " << tri[0] << " " << tri[1] << " " << tri[2] << endl;
    }
}

void writeSurfaceToVTK
(
    const word& name,
    const triSurf& surf,
    const List<DynList<scalar, 1> >& data
)
{
    OFstream file(name+".vtk");

    writeSurfaceToVTK(file, surf);

    //- write curvature fields
    const pointField& points = surf.points();
    file << "\n";
    file << "\nPOINT_DATA " << points.size() << "\n";

    file << "SCALARS Values double\n";
    file << "LOOKUP_TABLE default\n";
    forAll(points, pI)
    {
        file << data[pI][0] << " ";

        if( pI && pI % 5 == 0 )
            file << endl;
    }
}

void writeSurfaceToVTK
(
    const word& name,
    const triSurf& surf,
    const List<DynList<vector, 1> >& data
)
{
    OFstream file(name+".vtk");

    writeSurfaceToVTK(file, surf);

    //- write curvature fields
    const pointField& points = surf.points();
    file << "\n";
    file << "\nPOINT_DATA " << points.size() << "\n";

    file << "VECTORS Values double\n";
    file << "LOOKUP_TABLE default\n";
    forAll(points, pI)
    {
        const vector& v = data[pI][0];

        file << v[0] << " " << v[1] << " " << v[2] << " ";

        if( pI && pI % 5 == 0 )
            file << endl;
    }
}

void triSurfaceCurvatureEstimator::calculateEdgeCurvature()
{
    const pointField& points = surface_.points();
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& pointEdges = surface_.pointEdges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    edgePointCurvature_.setSize(points.size());
    boolList featureEdge(edges.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(edgePointCurvature_, i)
            edgePointCurvature_[i] = 0.0;

        //- mark feature edges
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(edgeFaces, eI)
        {
            if( edgeFaces.sizeOfRow(eI) != 2 )
            {
                featureEdge[eI] = false;
                continue;
            }

            if(
                 surface_[edgeFaces(eI, 0)].region() ==
                 surface_[edgeFaces(eI, 1)].region()
            )
            {
                featureEdge[eI] = false;
            }
            else
            {
                featureEdge[eI] = true;
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- loop through the points and calculate the curvature for points
        //- attached to two feature edges
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(pointEdges, pI)
        {
            DynList<label> features;
            forAllRow(pointEdges, pI, peI)
            {
                const label edgeI = pointEdges(pI, peI);
                if( featureEdge[edgeI] )
                    features.append(edgeI);
            }

            if( features.size() == 2 )
            {
                //- calculate the curvature and store it in the map
                vector e1 = edges[features[0]].vec(points);
                const scalar d1 = mag(e1);
                if( d1 > VSMALL )
                    e1 /= d1;
                vector e2 = edges[features[1]].vec(points);
                const scalar d2 = mag(e2);
                if( d2 > VSMALL )
                    e2 /= d2;

                scalar cs = e1 & e2;
                cs = Foam::min(1.0, cs);
                cs = Foam::max(1.0, cs);

                const scalar curv = Foam::acos(cs) / (0.5 * (d1+d2+VSMALL));
                edgePointCurvature_[pI] = Foam::mag(curv);
            }
        }
    }
}

void triSurfaceCurvatureEstimator::calculateSurfaceCurvatures()
{
    const VRWGraph& pointTriangles = surface_.pointFacets();

    const pointField& points = surface_.points();
    const VRWGraph& pointEdges = surface_.pointEdges();
    const edgeLongList& edges = surface_.edges();

    patchPositions_.setSize(surface_.size());
    gaussianCurvature_.setSize(points.size());
    meanCurvature_.setSize(points.size());
    maxCurvature_.setSize(points.size());
    minCurvature_.setSize(points.size());
    maxCurvatureVector_.setSize(points.size());
    minCurvatureVector_.setSize(points.size());

    List<DynList<label, 4> > pointPatches(points.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(pointTriangles, pointI)
    {
        std::map<label, DynList<label> > regionTriangles;
        Map<labelHashSet> otherLabels;
        Map<vector> normals;

        forAllRow(pointTriangles, pointI, ptI)
        {
            const label triI = pointTriangles(pointI, ptI);
            const label regionI = surface_[triI].region();

            if( !normals.found(regionI) )
                normals.insert(regionI, vector::zero);
            if( regionTriangles.find(regionI) != regionTriangles.end() )
                regionTriangles.insert
                (
                    std::make_pair
                    (
                        regionI,
                        DynList<label>()
                    )
                );
            if( !otherLabels.found(regionI) )
                otherLabels.insert(regionI, labelHashSet());

            regionTriangles[regionI].append(triI);

            forAll(surface_[triI], tpI)
            {
                const label pI = surface_[triI][tpI];

                if( pointI == pI )
                    continue;

                otherLabels[regionI].insert(pI);
                normals[regionI] += surface_[triI].normal(points);
            }
        }

        forAllConstIter(Map<labelHashSet>, otherLabels, it)
        {
            const labelHashSet& currLabels = it();

            labelHashSet additionalPoints;
            forAllConstIter(labelHashSet, currLabels, lit)
            {
                const label neiPointI = lit.key();
                const constRow pEdges = pointEdges[neiPointI];

                forAll(pEdges, peI)
                {
                    const label neiPointJ =
                        edges[pEdges[peI]].otherVertex(neiPointI);

                    if( neiPointJ == pointI )
                        continue;

                    bool correctPatch(false);
                    forAll(pointTriangles[neiPointJ], ntJ)
                    {
                        const label pJ = pointTriangles[neiPointJ][ntJ];
                        if( surface_[pJ].region() == it.key() )
                            correctPatch = true;
                    }

                    if( correctPatch )
                        additionalPoints.insert(neiPointJ);
                }
            }

            forAllConstIter(labelHashSet, additionalPoints, aIter)
                otherLabels[it.key()].insert(aIter.key());
        }

        forAllIter(Map<vector>, normals, nit)
            nit() /= (Foam::mag(nit()) + VSMALL);

        forAllConstIter(Map<labelHashSet>, otherLabels, it)
        {
            const labelHashSet& labels = it();

            const label regionI = it.key();

            //- store the patches
            pointPatches[pointI].append(regionI);

            //- find the position of point in the patch triangles
            const DynList<label>& rTriangles = regionTriangles[regionI];
            forAll(rTriangles, i)
            {
                const label tI = rTriangles[i];

                label pos(-1);
                forAll(surface_[tI], j)
                    if( surface_[tI][j] == pointI )
                    {
                        pos = j;
                        break;
                    }

                patchPositions_(tI, pos) = gaussianCurvature_[pointI].size();
            }

            //- store point coordinates
            DynList<point> op;

            forAllConstIter(labelHashSet, labels, lit)
                op.append(points[lit.key()]);

            //- fit the quadric patch to the surface
            quadricFitting qfit(points[pointI], normals[it.key()], op);

            # ifdef DEBUGCurvatureEstimator
            Info << "Point " << pointI << " in patch " << regionI
                << " has max curvature " << qfit.maxCurvature()
                << " and min curvature " << qfit.minCurvature() << endl;
            # endif

            //- store curvatures
            gaussianCurvature_[pointI].append(qfit.gaussianCurvature());
            meanCurvature_[pointI].append(qfit.meanCurvature());
            maxCurvature_[pointI].append(qfit.maxCurvature());
            minCurvature_[pointI].append(qfit.minCurvature());
            maxCurvatureVector_[pointI].append(qfit.maxCurvatureVector());
            minCurvatureVector_[pointI].append(qfit.minCurvatureVector());
        }
    }

    //- smooth curvatures using weighted Laplace
    List<DynList<scalar, 1> > smoothMinCurv(points.size());
    List<DynList<scalar, 1> > smoothMaxCurv(points.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        for(label iteration=0;iteration<2;++iteration)
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(pointEdges, pointI)
            {
                const constRow pEdges = pointEdges[pointI];

                //- find neighbouring points for each patch
                Map<DynList<label> > patchNeiPoints;
                forAll(pointPatches[pointI], ppI)
                    patchNeiPoints.insert
                    (
                        pointPatches[pointI][ppI],
                        DynList<label>()
                    );

                forAll(pEdges, peI)
                {
                    const edge& e = edges[pEdges[peI]];
                    const label neiPointI = e.otherVertex(pointI);

                    forAll(pointPatches[neiPointI], i)
                    {
                        const label patchI = pointPatches[neiPointI][i];

                        if( patchNeiPoints.found(patchI) )
                            patchNeiPoints[patchI].append(neiPointI);
                    }
                }

                smoothMinCurv[pointI].setSize(minCurvature_[pointI].size());
                smoothMaxCurv[pointI].setSize(maxCurvature_[pointI].size());

                //- update min curvature for all point patches
                forAll(minCurvature_[pointI], patchI)
                {
                    const label cPatch = pointPatches[pointI][patchI];

                    scalar minCurv(0.0);
                    scalar maxCurv(0.0);

                    const DynList<label>& neiPoints = patchNeiPoints[cPatch];

                    if( neiPoints.size() == 0 )
                    {
                        smoothMinCurv[pointI][patchI] =
                            minCurvature_[pointI][patchI];

                        smoothMaxCurv[pointI][patchI] =
                            maxCurvature_[pointI][patchI];
                    }

                    forAll(neiPoints, i)
                    {
                        const label neiPointI = neiPoints[i];

                        const label pos =
                            pointPatches[neiPointI].containsAtPosition(cPatch);

                        minCurv += minCurvature_[neiPointI][pos];
                        maxCurv += maxCurvature_[neiPointI][pos];
                    }

                    minCurv /= neiPoints.size();
                    maxCurv /= neiPoints.size();

                    //- store the value
                    smoothMinCurv[pointI][patchI] =
                        0.5 * (minCurv + minCurvature_[pointI][patchI]);

                    smoothMaxCurv[pointI][patchI] =
                        0.5 * (maxCurv + maxCurvature_[pointI][patchI]);
                }
            }

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp for schedule(static, 1)
            # endif
            forAll(minCurvature_, pointI)
            {
                forAll(minCurvature_[pointI], i)
                {
                    minCurvature_[pointI][i] = smoothMinCurv[pointI][i];
                    maxCurvature_[pointI][i] = smoothMaxCurv[pointI][i];
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- update Gaussian and mean curvatures
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(minCurvature_, pointI)
        {
            const DynList<scalar, 1>& minCurv = minCurvature_[pointI];
            const DynList<scalar, 1>& maxCurv = maxCurvature_[pointI];

            forAll(minCurv, i)
            {
                gaussianCurvature_[pointI][i] = minCurv[i] * maxCurv[i];
                meanCurvature_[pointI][i] = 0.5 * (minCurv[i] + maxCurv[i]);
            }
        }
    }

    # ifdef DEBUGCurvatureEstimator
    word name = "surfaceMeanCurv";
    writeSurfaceToVTK(name, surface_, meanCurvature_);
    writeSurfaceToVTK("surfaceGaussianCurv", surface_, gaussianCurvature_);
    writeSurfaceToVTK("surfaceMaxCurv", surface_, maxCurvature_);
    writeSurfaceToVTK("surfaceMinCurv", surface_, minCurvature_);
    writeSurfaceToVTK("surfaceMaxCurvVec", surface_, maxCurvatureVector_);
    writeSurfaceToVTK("surfaceMinCurvVec", surface_, minCurvatureVector_);
    # endif
}
/*
void triSurfaceCurvatureEstimator::calculateGaussianCurvature()
{
    gaussianCurvature_.setSize(surface_.patches().size());

    List<std::map<label, scalar> > magAreas(gaussianCurvature_.size());

    forAll(surface_, triI)
    {
        const labelledTri& tri = surface_[triI];
        const label region = tri.region();

        forAll(tri, pI)
        {
            gaussianCurvature_[region][tri[pI]] = 2 * M_PI;
            magAreas[region][tri[pI]] = 0.0;
        }
    }

    const pointField& points = surface_.points();
    const labelListList& pointTriangles = surface_.pointFaces();

    forAll(surface_, triI)
    {
        const labelledTri& tri = surface_[triI];
        const label region = tri.region();

        const scalar A = mag(tri.normal(points)) + VSMALL;

        const edgeList edges = tri.edges();
        vectorField ev(3);
        scalarField dv(3);
        forAll(edges, eI)
        {
            ev[eI] = edges[eI].vec(points);
            dv[eI] = mag(ev[eI]);
            if( dv[eI] > VSMALL )
                ev[eI] /= dv[eI];
        }

        forAll(tri, pI)
        {
            scalar cs = ev[(pI+2)%3] & ev[pI];
            cs = Foam::min(1.0, cs);
            cs = Foam::max(1.0, cs);

            gaussianCurvature_[region][tri[pI]] -= Foam::acos(-1.0 * cs);
            magAreas[region][tri[pI]] += A;
        }
    }

    //- calculate tge gaussian curvature for each vertex
    forAll(gaussianCurvature_, patchI)
    {
        std::map<label, scalar>& gc = gaussianCurvature_[patchI];
        std::map<label, scalar>& magSf = magAreas[patchI];

        for
        (
            std::map<label, scalar>::iterator it = gc.begin();
            it!=gc.end();
            ++it
        )
            it->second = 3.0 * it->second / magSf[it->first];
    }

    //- smooth the curvature variation
    const edgeList& edges = surface_.edges();
    const labelListList& pointEdges = surface_.pointEdges();
    forAll(meanCurvature_, patchI)
    {
        std::map<label, scalar>& gc = meanCurvature_[patchI];

        std::map<label, scalar> newValues;

        for(std::map<label, scalar>::iterator it=gc.begin();it!=gc.end();++it)
        {
            const labelList& pEdges = pointEdges[it->first];

            scalar maxCurv(-VSMALL), minCurv(VSMALL);
            label nNei(0);

            forAll(pEdges, peI)
            {
                const label ovI = edges[pEdges[peI]].otherVertex(it->first);

                if( gc.find(ovI) != gc.end() )
                {
                    maxCurv = Foam::max(maxCurv, gc[ovI]);
                    minCurv = Foam::min(minCurv, gc[ovI]);
                    ++nNei;
                }
            }

            if( nNei != 0 )
            {
                if( gc[it->first] > maxCurv )
                {
                    newValues.insert
                    (
                        std::make_pair(it->first, 0.5*(maxCurv+gc[it->first]))
                    );
                }
                else if( gc[it->first] >= minCurv )
                {
                    newValues.insert(*it);
                }
                else
                {
                    newValues.insert
                    (
                        std::make_pair(it->first, 0.5*(minCurv+gc[it->first]))
                    );
                }
            }
        }

        gc = newValues;
    }
}

void triSurfaceCurvatureEstimator::calculateMeanCurvature()
{
    const pointField& points = surface_.points();

    meanCurvature_.setSize(surface_.patches().size());

    List<Map<label> > nNeighbours(meanCurvature_.size());

    const edgeList& edges = surface_.edges();
    const labelListList& edgeFaces = surface_.edgeFaces();

    //- calculate edge angles
    forAll(edges, eI)
    {
        const labelList& eFaces = edgeFaces[eI];

        if( eFaces.size() != 2 )
            continue;
        if( surface_[eFaces[0]].region() != surface_[eFaces[1]].region() )
            continue;

        const label region = surface_[eFaces[0]].region();
        std::map<label, scalar>& mc = meanCurvature_[region];
        Map<label>& nn = nNeighbours[region];

        //- calculate normal vectors
        vector n0 = surface_[eFaces[0]].normal(points);
        const scalar A0 = mag(n0);
        if( A0 > VSMALL )
            n0 /= A0;

        vector n1 = surface_[eFaces[1]].normal(points);
        const scalar A1 = mag(n1);
        if( A1 > VSMALL )
            n1 /= A1;

        //- calculate the cosine of the angle between the two vectors
        scalar cs = n0 & n1;
        cs = Foam::min(1.0, cs);
        cs = Foam::max(-1.0, cs);

        //- calculate and normalize the edge vector
        vector ev = edges[eI].vec(points);
        const scalar length = mag(ev);
        if( length > VSMALL )
            ev /= length;

        //- calculate cross product of the normal vectors
        const vector t = n0 ^ n1;

        //- calculate the dot product between the edge vector and t
        scalar sn = t & ev;
        sn = Foam::min(1.0, sn);
        sn = Foam::max(-1.0, sn);

        //- calculate mean curvature over the edge
        scalar Hf(0.0);

        if( mag(sn) > VSMALL || mag(cs) > VSMALL )
        {
            //- calculate the angle
            const scalar angle = Foam::atan2(sn, cs);
            Hf += angle * length;
        }

        if( A0 > VSMALL || A1 > VSMALL )
        {
            Hf /= (A0 + A1 + VSMALL);
            Hf *= 3.0;
        }

        //- store the data
        const edge& e = edges[eI];
        if( mc.find(e.start()) == mc.end() )
        {
            mc.insert(std::make_pair(e.start(), 0.0));
            nn.insert(e.start(), 0);
        }
        mc[e.start()] += Hf;
        nn[e.start()] += 1;

        if( mc.find(e.end()) == mc.end() )
        {
            mc.insert(std::make_pair(e.end(), 0.0));
            nn.insert(e.end(), 0);
        }
        mc[e.end()] += Hf;
        nn[e.end()] += 1;
    }

    //- calculate the curvature for each vertex
    forAll(meanCurvature_, patchI)
    {
        std::map<label, scalar>& mc = meanCurvature_[patchI];
        Map<label>& nn = nNeighbours[patchI];

        forAllIter(Map<label>, nn, it)
        {
            if( it() == 0 )
            {
                mc[it.key()] = 0.0;
                continue;
            }

            mc[it.key()] = 0.5 * mc[it.key()] / it();
        }
    }

    //- smooth the curvature variation
    const labelListList& pointEdges = surface_.pointEdges();
    forAll(meanCurvature_, patchI)
    {
        std::map<label, scalar>& mc = meanCurvature_[patchI];

        std::map<label, scalar> newValues;

        for(std::map<label, scalar>::iterator it=mc.begin();it!=mc.end();++it)
        {
            const labelList& pEdges = pointEdges[it->first];

            scalar maxCurv(-VSMALL), minCurv(VSMALL);
            label nNei(0);

            forAll(pEdges, peI)
            {
                const label ovI = edges[pEdges[peI]].otherVertex(it->first);

                if( mc.find(ovI) != mc.end() )
                {
                    maxCurv = Foam::max(maxCurv, mc[ovI]);
                    minCurv = Foam::min(minCurv, mc[ovI]);
                    ++nNei;
                }
            }

            if( nNei != 0 )
            {
                if( mc[it->first] > maxCurv )
                {
                    newValues.insert
                    (
                        std::make_pair(it->first, 0.5*(maxCurv+mc[it->first]))
                    );
                }
                else if( mc[it->first] >= minCurv )
                {
                    newValues.insert(*it);
                }
                else
                {
                    newValues.insert
                    (
                        std::make_pair(it->first, 0.5*(minCurv+mc[it->first]))
                    );
                }
            }
        }

        mc = newValues;
    }
}

void triSurfaceCurvatureEstimator::calculateMinAndMaxCurvature()
{
    maxCurvature_.setSize(meanCurvature_.size());
    minCurvature_.setSize(meanCurvature_.size());

    forAll(maxCurvature_, patchI)
    {
        std::map<label, scalar>& mc = meanCurvature_[patchI];
        std::map<label, scalar>& gc = gaussianCurvature_[patchI];
        for
        (
            std::map<label, scalar>::const_iterator it=mc.begin();
            it!=mc.end();
            ++it
        )
        {
            const label pI = it->first;
            scalar disc = sqr(mc[pI]) - gc[pI];
            if( disc <  0 )
                disc = sqr(mc[pI]) + gc[pI];

            const scalar sqrtdisc = sqrt(disc);
            maxCurvature_[patchI][pI] = mc[pI] + sqrtdisc;
            minCurvature_[patchI][pI] = mc[pI] - sqrtdisc;
        }
    }
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
