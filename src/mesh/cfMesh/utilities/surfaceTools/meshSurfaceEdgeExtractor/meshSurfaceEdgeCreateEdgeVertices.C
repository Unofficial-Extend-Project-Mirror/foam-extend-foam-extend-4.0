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
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEdgeExtractor.H"
#include "meshOctree.H"
#include "Map.H"

#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEdgeExtractor::createEdgeVertices()
{
    pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    const edgeList& edges = mesh_.addressingData().edges();
    const VRWGraph& faceEdges = mesh_.addressingData().faceEdges();

    Map<label> newEdgePoint;

    const label nIntFaces = mesh_.nInternalFaces();
    const label nFaces = faces.size();

    nPoints_ = points.size();

    for(label faceI=nIntFaces;faceI<nFaces;++faceI)
    {
        const face& f = faces[faceI];

        forAll(f, pI)
        {
            const label edgeI = faceEdges(faceI, pI);

            if( newEdgePoint.found(edgeI) ) continue;

            const label s = f[pI];
            const label e = f.nextLabel(pI);

            if( !pointRegions_.sizeOfRow(s) || !pointRegions_.sizeOfRow(e) )
            {
                Warning << "Boundary vertices " << s << " and " << e
                    << " are not mapped to the boundary!" << endl;

                continue;
            }

            if( pointRegions_(s, 0) != pointRegions_(e, 0) )
            {
                point newP;
                scalar distSq;
                label nse;

                FixedList<point, 2> edgePoints;
                FixedList<label, 2> patches;

                edgePoints[0] = points[s];
                edgePoints[1] = points[e];
                patches[0] = pointRegions_(s, 0);
                patches[1] = pointRegions_(e, 0);

                const bool found =
                meshOctree_.findNearestPointToEdge
                (
                    newP,
                    distSq,
                    nse,
                    edgePoints,
                    patches
                );

                if( found )
                {
                    points.append(newP);
                }
                else
                {
                    points.append
                    (
                        edges[faceEdges(faceI, pI)].centre(points)
                    );
                }

                pointRegions_.appendList(patches);

                newEdgePoint.insert(edgeI, nPoints_);
                ++nPoints_;
            }
        }
    }

    points.setSize(nPoints_);

    //- create new faces
    DynList<label> newF;
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        newF.clear();

        forAll(f, eI)
        {
            newF.append(f[eI]);
            if( newEdgePoint.found(faceEdges(faceI, eI)) )
                newF.append(newEdgePoint[faceEdges(faceI, eI)]);
        }

        if( newF.size() > f.size() )
        {
            //- face must be changed
            face& mf = const_cast<face&>(f);
            mf.setSize(newF.size());
            forAll(mf, pI)
                mf[pI] = newF[pI];
        }
    }

    mesh_.clearAddressingData();

    Info << "Finished creating mesh edges" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
