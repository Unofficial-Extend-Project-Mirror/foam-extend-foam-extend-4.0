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

#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"
#include "labelLongList.H"
#include "boolList.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfacePartitioner::calculatePatchAddressing()
{
    calculateCornersAndAddressing();

    calculatePatchPatches();

    calculateEdgeGroups();

    calculatePatchToEdgeGroups();

    calculateEdgeGroupsToCorners();
}

void triSurfacePartitioner::calculateCornersAndAddressing()
{
    const VRWGraph& pointFaces = surface_.pointFacets();
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    //- find the number of feature edges connected to each surface node
    labelList nEdgesAtNode(surface_.points().size(), 0);
    forAll(edgeFaces, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 )
            continue;

        const label sPatch = surface_[edgeFaces(eI, 0)].region();
        const label ePatch = surface_[edgeFaces(eI, 1)].region();

        if( sPatch != ePatch )
        {
            const edge& e = edges[eI];
            ++nEdgesAtNode[e.start()];
            ++nEdgesAtNode[e.end()];
        }
    }

    //- count the number of feature edges connected to each surface point
    //- corners must have 3 or more edges attached to them
    label nCorners(0);
    forAll(nEdgesAtNode, pI)
    {
        if( nEdgesAtNode[pI] < 3 )
            continue;

        ++nCorners;
    }

    corners_.setSize(nCorners);
    cornerPatches_.setSize(nCorners);
    nCorners = 0;

    //- store corner data
    DynList<label> patches;
    forAll(pointFaces, pointI)
    {
        if( nEdgesAtNode[pointI] < 3 )
            continue;

        patches.clear();
        forAllRow(pointFaces, pointI, pfI)
            patches.appendIfNotIn(surface_[pointFaces(pointI, pfI)].region());

        corners_[nCorners] = pointI;
        cornerPatches_[nCorners] = patches;
        ++nCorners;
    }
}

void triSurfacePartitioner::calculatePatchPatches()
{
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    forAll(edgeFaces, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 )
        {
            Warning << "Surface is not a manifold!!" << endl;
            continue;
        }

        const label sPatch = surface_[edgeFaces(eI, 0)].region();
        const label ePatch = surface_[edgeFaces(eI, 1)].region();

        if( sPatch != ePatch )
        {
            patchPatches_[sPatch].insert(ePatch);
            patchPatches_[ePatch].insert(sPatch);
        }
    }
}

void triSurfacePartitioner::calculateEdgeGroups()
{
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();
    const VRWGraph& pointEdges = surface_.pointEdges();


    //- make all feature edges
    boolList featureEdge(edgeFaces.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(edgeFaces, eI)
    {
        DynList<label> patches;
        forAllRow(edgeFaces, eI, efI)
            patches.appendIfNotIn(surface_[edgeFaces(eI, efI)].region());

        if( patches.size() > 1 )
            featureEdge[eI] = true;
    }

    //- create a set containing corners for fast searching
    labelHashSet corners;
    forAll(corners_, i)
        corners.insert(corners_[i]);

    edgeGroups_.setSize(edgeFaces.size());
    edgeGroups_ = -1;

    label nGroups(0);
    forAll(featureEdge, eI)
    {
        if( !featureEdge[eI] )
            continue;
        if( edgeGroups_[eI] != -1 )
            continue;

        labelLongList front;
        front.append(eI);
        edgeGroups_[eI] = nGroups;
        featureEdge[eI] = false;

        while( front.size() )
        {
            const label eLabel = front.removeLastElement();
            const edge& e = edges[eLabel];

            for(label pI=0;pI<2;++pI)
            {
                const label pointI = e[pI];

                if( corners.found(pointI) )
                    continue;

                forAllRow(pointEdges, pointI, peI)
                {
                    const label eJ = pointEdges(pointI, peI);

                    if( featureEdge[eJ] && (edgeGroups_[eJ] == -1) )
                    {
                        edgeGroups_[eJ] = nGroups;
                        featureEdge[eJ] = false;
                        front.append(eJ);
                    }
                }
            }
        }

        ++nGroups;
    }

    Info << nGroups << " edge groups found!" << endl;

    edgeGroupEdgeGroups_.clear();
    edgeGroupEdgeGroups_.setSize(nGroups);
}

void triSurfacePartitioner::calculatePatchToEdgeGroups()
{
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    forAll(edgeFaces, eI)
    {
        if( edgeGroups_[eI] < 0 )
            continue;

        DynList<label> patches;
        forAllRow(edgeFaces, eI, efI)
            patches.appendIfNotIn(surface_[edgeFaces(eI, efI)].region());

        forAll(patches, i)
        {
            const label patchI = patches[i];

            for(label j=i+1;j<patches.size();++j)
            {
                const label patchJ = patches[j];

                const std::pair<label, label> pp
                (
                    Foam::min(patchI, patchJ),
                    Foam::max(patchI, patchJ)
                );

                patchesEdgeGroups_[pp].insert(edgeGroups_[eI]);
            }
        }
    }
}

void triSurfacePartitioner::calculateEdgeGroupsToCorners()
{
    const VRWGraph& pointEdges = surface_.pointEdges();

    forAll(corners_, cornerI)
    {
        DynList<label> edgeGroupsAtCorner;
        const label pointI = corners_[cornerI];

        forAllRow(pointEdges, pointI, peI)
            edgeGroupsAtCorner.appendIfNotIn
            (
                edgeGroups_[pointEdges(pointI, peI)]
            );

        forAll(edgeGroupsAtCorner, i)
        {
            const label epI = edgeGroupsAtCorner[i];
            if( epI < 0 )
                continue;
            for(label j=i+1;j<edgeGroupsAtCorner.size();++j)
            {
                const label epJ = edgeGroupsAtCorner[j];
                if( epJ < 0 )
                    continue;

                std::pair<label, label> ep
                (
                    Foam::min(epI, epJ),
                    Foam::max(epI, epJ)
                );

                //- create edgepartition - edge partitions addressing
                edgeGroupEdgeGroups_[ep.first].insert(ep.second);
                edgeGroupEdgeGroups_[ep.second].insert(ep.first);

                edgeGroupsCorners_[ep].insert(cornerI);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
