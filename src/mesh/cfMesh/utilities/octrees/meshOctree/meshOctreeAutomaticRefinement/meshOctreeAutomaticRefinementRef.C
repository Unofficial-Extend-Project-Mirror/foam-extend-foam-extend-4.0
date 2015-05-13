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

#include "meshOctreeAutomaticRefinement.H"
#include "triSurface.H"
#include "demandDrivenData.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceCurvatureEstimator.H"
#include "meshOctreeAddressing.H"
#include "triSurf.H"
#include "helperFunctions.H"
#include "meshOctreeModifier.H"

#include "Map.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGAutoRef

# ifdef DEBUGAutoRef
#include "pointSet.H"
#include "IOdictionary.H"
#include "objectRegistry.H"
#include "Time.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeAutomaticRefinement::activateHexRefinement()
{
    hexRefinement_ = true;
}

void meshOctreeAutomaticRefinement::automaticRefinement()
{
    Info << "Performing automatic refinement" << endl;

    if( !maxRefLevel_ )
        return;

    curvatureRefinement();

    proximityRefinement();

    Info << "Finished with automatic refinement" << endl;
}

bool meshOctreeAutomaticRefinement::curvatureRefinement()
{
    List<direction> refineBox(octree_.numberOfLeaves(), direction(0));
    labelLongList refinementCandidates;
    forAll(refineBox, i)
        refinementCandidates.append(i);
    while( refineBasedOnCurvature(refineBox, refinementCandidates) )
    {
        refineSelectedBoxes(refineBox, refinementCandidates);
    }

    return false;
}

bool meshOctreeAutomaticRefinement::proximityRefinement()
{
    bool refine(false);
    List<direction> refineBox(octree_.numberOfLeaves(), direction(0));
    labelLongList refinementCandidates;
    forAll(refineBox, i)
        refinementCandidates.append(i);
    while( refineBasedOnContainedCorners(refineBox, refinementCandidates) )
    {
        refineSelectedBoxes(refineBox, refinementCandidates);
        refine = true;
    }

    refinementCandidates.clear();
    forAll(refineBox, i)
        refinementCandidates.append(i);
    while( refineBasedOnContainedPartitions(refineBox, refinementCandidates) )
    {
        refineSelectedBoxes(refineBox, refinementCandidates);
        refine = true;
    }

    refinementCandidates.clear();
    forAll(refineBox, i)
        refinementCandidates.append(i);
    while( refineBasedOnProximityTests(refineBox, refinementCandidates) )
    {
        refineSelectedBoxes(refineBox, refinementCandidates);
        refine = true;
    }

    return refine;
}

bool meshOctreeAutomaticRefinement::refineBasedOnContainedCorners
(
    List<direction>& refineBox,
    const labelLongList& refCandidates
)
{
    meshOctreeModifier octreeModifier(octree_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();
    const boundBox& rootBox = octree_.rootBox();
    const triSurf& surface = octree_.surface();
    const pointField& points = surface.points();
    const triSurfacePartitioner& sPart = this->partitioner();

    //- find leaves which contains corner nodes
    labelList cornerInLeaf(refineBox.size(), -1);
    const labelList& corners = sPart.corners();

    label nMarked(0);

    forAll(corners, cornerI)
    {
        const label cLabel =
        octree_.findLeafContainingVertex(points[corners[cornerI]]);

        if( cLabel < 0 )
            continue;

        if( cornerInLeaf[cLabel] != -1 )
        {
            // refine this box because it already contains some other corner
            ++nMarked;
            refineBox[cLabel] = 1;
        }
        else
        {
            cornerInLeaf[cLabel] = corners[cornerI];
        }
    }

    DynList<label> leavesInBox;
    # ifdef USE_OMP
    # pragma omp parallel for if( refCandidates.size() > 1000 ) \
    private(leavesInBox) shared(cornerInLeaf) \
    reduction(+ : nMarked) schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];
        if( leaves[leafI]->level() >= maxRefLevel_ )
            continue;
        if( cornerInLeaf[leafI] == -1 )
            continue;

        // check if there exist some corners in the neighbour boxes
        // refine the box if some corners are found in the neighbouring boxes
        const point c = leaves[leafI]->centre(rootBox);
        const scalar r = 1.732 * leaves[leafI]->size(rootBox);

        boundBox bb(c - point(r, r, r), c + point(r, r, r));

        leavesInBox.clear();
        octree_.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label nei = leavesInBox[i];

            if( nei < 0 )
                continue;

            if( nei == leafI )
                continue;
            if( cornerInLeaf[nei] == -1 )
                continue;

            if( mag(points[cornerInLeaf[nei]] - c) < r )
            {
                ++nMarked;
                refineBox[nei] = 1;
                refineBox[leafI] = 1;
                break;
            }
        }
    }

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxes marked by the corner criteria" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnContainedPartitions
(
    List<direction>& refineBox,
    const labelLongList& refCandidates
)
{
    const boundBox& rootBox = octree_.rootBox();
    const triSurfacePartitioner& sPart = this->partitioner();

    //- find leaves which contains corner nodes
    const List<labelHashSet>& pPatches = sPart.patchPatches();
    const labelList& edgeGroups = sPart.edgeGroups();
    const List<labelHashSet>& eNeiGroups = sPart.edgeGroupEdgeGroups();

    # ifdef DEBUGAutoRef
    Info << "pPart " << pPart << endl;
    # endif

    label nMarked(0);

    meshOctreeModifier octreeModifier(octree_);
    const triSurf& surf = octree_.surface();
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    DynList<label> patches, eGroups, helper;
    # ifdef USE_OMP
    # pragma omp parallel for if( refCandidates.size() > 1000 ) \
    private(patches, eGroups, helper) \
    reduction(+ : nMarked) schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];
        if( !leaves[leafI]->hasContainedElements() )
            continue;
        if( leaves[leafI]->level() >= maxRefLevel_ )
            continue;

        const meshOctreeCubeBasic& oc = *leaves[leafI];
        const point c = oc.centre(rootBox);
        const scalar s = 1.733 * oc.size(rootBox);
        const boundBox bb(c - point(s, s, s), c + point(s, s, s));

        //- find triangle patches contained in this box
        octree_.findTrianglesInBox(bb, helper);
        patches.clear();
        forAll(helper, i)
            patches.appendIfNotIn(surf[helper[i]].region());

        //- find edge partitions contained in this box
        helper.clear();
        octree_.findEdgesInBox(bb, helper);
        eGroups.clear();
        forAll(helper, i)
            eGroups.appendIfNotIn(edgeGroups[helper[i]]);

        # ifdef DEBUGAutoRef
        Info << "patches for leaf " << leafI << " are " << patches << endl;
        # endif

        bool refine(false);
        forAll(patches, patchI)
        {
            for(label patchJ=(patchI+1);patchJ<patches.size();++patchJ)
                if( !pPatches[patches[patchI]].found(patches[patchJ]) )
                {
                    # ifdef DEBUGAutoRef
                    Info << "2.Here" << endl;
                    # endif

                    refine = true;
                    break;
                }
        }

        forAll(eGroups, egI)
        {
            for(label egJ=egI+1;egJ<eGroups.size();++egJ)
                if( !eNeiGroups[eGroups[egI]].found(eGroups[egJ]) )
                {
                    refine = true;
                    break;
                }
        }

        if( refine )
        {
            # ifdef DEBUGAutoRef
            Info << "Selecting leaf " << leafI
                << " for auto-refinement" << endl;
            # endif

            ++nMarked;
            refineBox[leafI] = 1;
        }
    }

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxed marked by partitioning criteria" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnCurvature
(
    List<direction>& refineBox,
    const labelLongList& refCandidates
)
{
    const triSurfaceCurvatureEstimator& curv = curvature();

    const boundBox& rootBox = octree_.rootBox();

    label nMarked(0);
    DynList<label> containedTrias;
    # ifdef USE_OMP
    # pragma omp parallel for if( refCandidates.size() > 10000 ) \
    private(containedTrias) \
    reduction(+ : nMarked) schedule(dynamic, 100)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel_ )
            continue;

        //- search for the minimum curvature radius at surface triangles
        octree_.containedTriangles(leafI, containedTrias);

        scalar maxCurv(0.0);
        forAll(containedTrias, i)
            maxCurv =
                Foam::max
                (
                    maxCurv,
                    mag(curv.meanCurvatureAtTriangle(containedTrias[i]))
                );

        //- check the edge curvature
        if( octree_.hasContainedEdges(leafI) )
        {
            octree_.containedEdges(leafI, containedTrias);

            forAll(containedTrias, i)
            {
                maxCurv =
                    Foam::max
                    (
                        maxCurv,
                        mag(curv.curvatureAtEdge(containedTrias[i]))
                    );
            }
        }

        if( oc.size(rootBox) > 0.2835 / (maxCurv + SMALL) )
        {
            refineBox[leafI] = 1;
            ++nMarked;
        }
    }

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxes marked by curvature criteria!" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnProximityTests
(
    List<direction>& refineBox,
    const labelLongList& refCandidates
)
{
    const boundBox& rootBox = octree_.rootBox();
    const triSurf& surf = octree_.surface();

    label nMarked(0);
    DynList<label> helper;
    # ifdef USE_OMP
    # pragma omp parallel for if( refCandidates.size() > 1000 ) \
    private(helper) reduction(+ : nMarked) schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel_ )
            continue;

        const point c = oc.centre(rootBox);
        const scalar s = 1.732 * oc.size(rootBox);
        boundBox bb(c - point(s, s, s), c + point(s, s, s));

        labelHashSet triaInRange(100), edgesInRange(100);
        //- find triangles in range
        helper.clear();
        octree_.findTrianglesInBox(bb, helper);
        forAll(helper, i)
            triaInRange.insert(helper[i]);

        //- find edges contained in the neighbourhood
        helper.clear();
        octree_.findEdgesInBox(bb, helper);
        forAll(helper, i)
            edgesInRange.insert(helper[i]);

        //- refine boxes with more than two face groups
        if
        (
            (help::numberOfFaceGroups(triaInRange, c, s, surf) > 1) ||
            (help::numberOfEdgeGroups(edgesInRange, c, s, surf) > 1)
        )
        {
            ++nMarked;
            refineBox[leafI] = 1;
        }
    }

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxed marked by proximity criteria" << endl;

    if( nMarked != 0 )
        return true;

    return false;
}

void meshOctreeAutomaticRefinement::refineSelectedBoxes
(
    List<direction>& refineBox,
    labelLongList& refCandidates
)
{
    deleteDemandDrivenData(octreeAddressingPtr_);

    meshOctreeModifier octreeModifier(octree_);
    LongList<meshOctreeCube*> leaves = octreeModifier.leavesAccess();

    octreeModifier.markAdditionalLayers(refineBox, 1);
    octreeModifier.refineSelectedBoxes(refineBox, hexRefinement_);

    //- find the cubes which have been marked for refinement
    LongList<meshOctreeCubeCoordinates> refinedCubes;
    forAll(refineBox, i)
    {
        if( refineBox[i] )
            refinedCubes.append(leaves[i]->coordinates());
    }
    leaves.setSize(0);

    //- perform load distribution in case od parallel runs
    octreeModifier.loadDistribution();

    //- communicate the cubes selected for refinement with other processors
    LongList<meshOctreeCubeCoordinates> receivedCoordinates;
    if( Pstream::parRun() )
    {
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            refinedCubes,
            receivedCoordinates
        );
    }

    forAll(refinedCubes, i)
        receivedCoordinates.append(refinedCubes[i]);
    refinedCubes.setSize(0);

    //- find the cubes which shall checked in the next iteration
    refCandidates.clear();
    forAll(receivedCoordinates, i)
    {
        const meshOctreeCubeCoordinates& cc = receivedCoordinates[i];

        for(label scI=0;scI<8;++scI)
        {
            const meshOctreeCubeCoordinates child = cc.refineForPosition(scI);

            meshOctreeCube* oc = octreeModifier.findCubeForPosition(child);

            if( !oc || !oc->isLeaf() )
                continue;

            refCandidates.append(oc->cubeLabel());
        }
    }

    refineBox.setSize(octree_.numberOfLeaves());
    refineBox = direction(0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
