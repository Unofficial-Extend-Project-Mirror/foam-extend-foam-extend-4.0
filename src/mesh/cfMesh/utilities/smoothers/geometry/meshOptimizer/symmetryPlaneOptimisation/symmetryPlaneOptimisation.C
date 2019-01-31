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
#include "symmetryPlaneOptimisation.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"
#include "polyMeshGenChecks.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void symmetryPlaneOptimisation::detectSymmetryPlanes()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    symmetryPlanes_.clear();

    typedef std::map<label, std::pair<vector, label> > mapType;
    mapType centreSum, normalSum;

    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].patchType() == "symmetryPlane" )
        {
            std::pair<vector, label>& cs = centreSum[patchI];
            cs = std::pair<vector, label>(vector::zero, 0);

            std::pair<vector, label>& ns = normalSum[patchI];
            ns = std::pair<vector, label>(vector::zero, 0);

            const label start = boundaries[patchI].patchStart();
            const label end = start + boundaries[patchI].patchSize();
            for(label faceI=start;faceI<end;++faceI)
            {
                cs.first += faces[faceI].centre(points);
                ns.first += faces[faceI].normal(points);
            }

            cs.second = ns.second = boundaries[patchI].patchSize();
        }
    }

    if( Pstream::parRun() )
    {
        //- sum up all normals and centres of all processors
        //- every symmetry plane patch must be present on all processors
        forAllIter(mapType, centreSum, pIter)
        {
            std::pair<vector, label>& cs = pIter->second;
            reduce(cs.second, sumOp<label>());
            reduce(cs.first, sumOp<vector>());

            std::pair<vector, label>& ns = normalSum[pIter->first];
            reduce(ns.first, sumOp<vector>());
            ns.second = cs.second;
        }
    }

    //- create planes corresponding to each symmetry plane
    forAllConstIter(mapType, centreSum, it)
    {
        const point c = it->second.first / it->second.second;

        const std::pair<vector, label>& ns = normalSum[it->first];
        const vector n = ns.first / ns.second;

        symmetryPlanes_.insert(std::make_pair(it->first, plane(c, n)));
    }
}

bool symmetryPlaneOptimisation::pointInPlanes(VRWGraph& pointInPlanes) const
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    pointInPlanes.clear();
    pointInPlanes.setSize(points.size());

    bool foundProblematic(false);

    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].patchType() == "symmetryPlane" )
        {
            const label start = boundaries[patchI].patchStart();
            const label end = start + boundaries[patchI].patchSize();
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                const point c = f.centre(points);
                scalar maxDist(0.0);
                forAll(f, pI)
                    maxDist = max(maxDist, mag(points[f[pI]] - c));

                forAll(f, pI)
                {
                    std::map<label, plane>::const_iterator it =
                        symmetryPlanes_.find(patchI);
                    if( it != symmetryPlanes_.end() )
                    {
                        const scalar d = it->second.distance(points[f[pI]]);

                        if( d > 0.5 * maxDist )
                            foundProblematic = true;
                    }

                    pointInPlanes.appendIfNotIn(f[pI], patchI);
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalPointAddressing();
        const VRWGraph& pointAtProcs = mesh_.addressingData().pointAtProcs();
        const DynList<label>& neiProcs =
            mesh_.addressingData().pointNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            if( pointInPlanes.sizeOfRow(pointI) == 0 )
                continue;

            forAllRow(pointAtProcs, pointI, i)
            {
                const label neiProc = pointAtProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& dataToSend = exchangeData[neiProc];

                dataToSend.append(it.key());
                dataToSend.append(pointInPlanes.sizeOfRow(pointI));
                forAllRow(pointInPlanes, pointI, pipI)
                    dataToSend.append(pointInPlanes(pointI, pipI));
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        for(label counter=0;counter<receivedData.size();)
        {
            const label pointI = globalToLocal[receivedData[counter++]];

            const label size = receivedData[counter++];
            for(label i=0;i<size;++i)
                pointInPlanes.appendIfNotIn(pointI, receivedData[counter++]);
        }
    }

    reduce(foundProblematic, maxOp<bool>());

    return foundProblematic;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
symmetryPlaneOptimisation::symmetryPlaneOptimisation(polyMeshGen& mesh)
:
    mesh_(mesh)
{
    detectSymmetryPlanes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

symmetryPlaneOptimisation::~symmetryPlaneOptimisation()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void symmetryPlaneOptimisation::optimizeSymmetryPlanes()
{
    pointFieldPMG& points = mesh_.points();

    VRWGraph pointInPlane;
    if( pointInPlanes(pointInPlane) )
    {
        Warning
            << "Detected large deviation from some symmetry planes."
            << "Please check your settings and ensure that all your patches"
            << " with type symmetryPlane are planar. If there exists a patch"
            << " that is not planar please use the symmetry type, instead."
            << " Skipping symmetryPlane optimisation..." << endl;

        return;
    }

    forAll(pointInPlane, pointI)
    {
        const label nPlanes = pointInPlane.sizeOfRow(pointI);

        if( nPlanes > 3 )
        {
            WarningIn
            (
                "void symmetryPlaneOptimisation::optimizeSymmetryPlanes()"
            ) << "Point " << pointI << " is in more than three symmetry"
              << " planes. Cannot move it" << endl;

            continue;
        }

        point& p = points[pointI];

        if( nPlanes == 1 )
        {
            //- point is in a plane
            std::map<label, plane>::const_iterator it =
                symmetryPlanes_.find(pointInPlane(pointI, 0));

            p = it->second.nearestPoint(p);
        }
        else if( nPlanes == 2 )
        {
            //- point is at the edge between two planes
            const plane& pl0 =
                symmetryPlanes_.find(pointInPlane(pointI, 0))->second;
            const plane& pl1 =
                symmetryPlanes_.find(pointInPlane(pointI, 1))->second;

            const plane::ray il = pl0.planeIntersect(pl1);

            vector n = il.dir();
            n /= (mag(n) + VSMALL);

            const scalar lProj = (p - il.refPoint()) & n;

            p = il.refPoint() + lProj * n;
        }
        else if( nPlanes == 3 )
        {
            //- points is a corner between three planes
            const plane& pl0 =
                symmetryPlanes_.find(pointInPlane(pointI, 0))->second;
            const plane& pl1 =
                symmetryPlanes_.find(pointInPlane(pointI, 1))->second;
            const plane& pl2 =
                symmetryPlanes_.find(pointInPlane(pointI, 2))->second;

            p = pl0.planePlaneIntersect(pl1, pl2);
        }
    }

    labelHashSet badFaces;
    polyMeshGenChecks::checkFacePyramids(mesh_, false, VSMALL, &badFaces);

    if( badFaces.size() )
    {
        WarningIn
        (
            "void symmetryPlaneOptimisation::optimizeSymmetryPlanes()"
        ) << "Bad quality or inverted faces found in the mesh" << endl;

        const label badFacesId = mesh_.addFaceSubset("invalidFaces");
        forAllConstIter(labelHashSet, badFaces, it)
            mesh_.addFaceToSubset(badFacesId, it.key());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
