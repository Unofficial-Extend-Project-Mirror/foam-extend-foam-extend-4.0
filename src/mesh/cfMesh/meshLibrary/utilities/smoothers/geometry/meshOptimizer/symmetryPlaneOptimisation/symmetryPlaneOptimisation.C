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
#include "meshOptimizer.H"

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
        const point n = ns.first / ns.second;

        symmetryPlanes_.insert(std::make_pair(it->first, plane(c, n)));
    }
}

void symmetryPlaneOptimisation::pointInPlanes(VRWGraph& pointInPlanes) const
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    pointInPlanes.clear();
    pointInPlanes.setSize(points.size());

    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].patchType() == "symmetryPlane" )
        {
            const label start = boundaries[patchI].patchStart();
            const label end = start + boundaries[patchI].patchSize();
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                forAll(f, pI)
                    pointInPlanes.appendIfNotIn(f[pI], patchI);
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
    pointInPlanes(pointInPlane);

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
        }

        point& p = points[pointI];
        vector disp(vector::zero);
        for(label plI=0;plI<nPlanes;++plI)
        {
            //- point is in a plane
            std::map<label, plane>::const_iterator it =
                symmetryPlanes_.find(pointInPlane(pointI, 0));

            const point newP = it->second.nearestPoint(points[pointI]);
            disp += newP - p;
        }

        p += disp;
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
