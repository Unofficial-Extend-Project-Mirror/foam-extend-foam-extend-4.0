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

#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "boundBox.H"
#include "demandDrivenData.H"
#include "objectRefinementList.H"
#include "VRWGraph.H"
#include "meshOctreeModifier.H"
#include "helperFunctions.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define OCTREETiming
//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCreator::refineBoundary()
{
    meshOctreeModifier octreeModifier(octree_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    //- refine DATA boxes to the given level
    Info << "Refining boundary boxes to the given size" << endl;

    label nMarked;
    do
    {
        nMarked = 0;

        # ifdef OCTREETiming
        const scalar startIter = omp_get_wtime();
        # endif

        List<direction> refineCubes(leaves.size(), direction(0));
        labelList nLayers(leaves.size(), 0);
        List<direction> targetLevel(leaves.size(), direction(0));
        bool useNLayers(false);

        //- select boxes which need to be refined
        # ifdef USE_OMP
        # pragma omp parallel for reduction(+ : nMarked) \
        schedule(dynamic, Foam::min(20, leaves.size()/omp_get_num_threads()+1))
        # endif
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];
            # ifdef DEBUGSearch
            Info << "Checking leaf " << oc << endl;
            Info << "Leaf has elements " << oc.hasContainedElements() << endl;
            # endif

            if( oc.hasContainedElements() )
            {
                const label elRowI = oc.containedElements();
                const VRWGraph& containedTriangles =
                    oc.slotPtr()->containedTriangles_;

                const scalar cs = oc.size(octree_.rootBox());

                bool refine(false);
                forAllRow(containedTriangles, elRowI, tI)
                {
                    const label triI = containedTriangles(elRowI, tI);

                    if( surfRefLevel_[triI] > oc.level() )
                    {
                        refine = true;
                    }

                    if( surfRefThickness_[triI] > VSMALL )
                    {
                        useNLayers = true;

                        nLayers[leafI] =
                            Foam::max
                            (
                                nLayers[leafI],
                                Foam::max(label(surfRefThickness_[triI]/cs), 1)
                            );

                        targetLevel[leafI] =
                            Foam::max(targetLevel[leafI], surfRefLevel_[triI]);
                    }
                }

                if( refine )
                {
                    refineCubes[leafI] = 1;
                    ++nMarked;
                }
            }
        }

        //- mark additional boxes for refinement to achieve
        //- correct refinement distance
        reduce(useNLayers, maxOp<label>());
        if( useNLayers )
        {
            nMarked +=
                octreeModifier.markAdditionalLayers
                (
                    refineCubes,
                    nLayers,
                    targetLevel
                );
        }

        //- refine boxes
        octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);

        # ifdef OCTREETiming
        const scalar refTime = omp_get_wtime();
        Info << "Time for refinement " << (refTime-startIter) << endl;
        # endif

        if( Pstream::parRun() )
        {
            reduce(nMarked, sumOp<label>());
            if( nMarked )
            {
                octreeModifier.distributeLeavesToProcessors();

                # ifdef OCTREETiming
                const scalar distTime = omp_get_wtime();
                Info << "Time for distributing data to processors "
                     << (distTime-refTime) << endl;
                # endif

                loadDistribution();

                # ifdef OCTREETiming
                Info << "Time for load distribution "
                    << (omp_get_wtime()-distTime) << endl;
                # endif
            }
        }

    } while( nMarked );

    Info << "Finished refining boundary boxes" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCreator::refineBoxesContainedInObjects()
{
    if( !meshDictPtr_ || !meshDictPtr_->found("objectRefinements") )
    {
        return;
    }

    Info << "Refining boxes inside objects" << endl;
    objectRefinementList refObjects;

    // Read objects
    if( meshDictPtr_->isDict("objectRefinements") )
    {
        const dictionary& dict = meshDictPtr_->subDict("objectRefinements");
        const wordList objectNames = dict.toc();

        refObjects.setSize(objectNames.size());

        forAll(refObjects, objectI)
        {
            const entry& objectEntry =
                dict.lookupEntry(objectNames[objectI], false, false);

            refObjects.set
            (
                objectI,
                objectRefinement::New
                (
                    objectEntry.keyword(),
                    objectEntry.dict()
                )
            );
        }
    }
    else
    {
        Istream& is = meshDictPtr_->lookup("objectRefinements");

        PtrList<entry> objectEntries(is);
        refObjects.setSize(objectEntries.size());

        forAll(refObjects, objectI)
        {
            refObjects.set
            (
                objectI,
                objectRefinement::New
                (
                    objectEntries[objectI].keyword(),
                    objectEntries[objectI].dict()
                )
            );
        }

        objectEntries.clear();
    }

    scalar s(readScalar(meshDictPtr_->lookup("maxCellSize")));

    List<direction> refLevels(refObjects.size(), globalRefLevel_);
    scalarList refThickness(refObjects.size(), 0.0);

    forAll(refThickness, oI)
        refThickness[oI] = refObjects[oI].refinementThickness();

    label nMarked;
    do
    {
        nMarked = 0;
        forAll(refObjects, oI)
        {
            if( refObjects[oI].cellSize() <= s * (1.+SMALL) )
            {
                ++nMarked;
                ++refLevels[oI];
            }
        }

        s /= 2.0;

    } while( nMarked != 0 );

    forAll(refLevels, i)
        Info << "Ref level for object " << refObjects[i].name()
            << " is " << label(refLevels[i]) << endl;

    if( octree_.neiProcs().size() )
        forAll(refObjects, oI)
        {
            label l = refLevels[oI];
            reduce(l, maxOp<label>());
            refLevels[oI] = l;
        }

    //- start refining boxes inside the objects
    const boundBox& rootBox = octree_.rootBox();
    meshOctreeModifier octreeModifier(octree_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    do
    {
        # ifdef OCTREETiming
        const scalar startIter = omp_get_wtime();
        # endif

        nMarked = 0;

        List<direction> refineCubes(leaves.size(), direction(0));
        labelList nLayers(leaves.size(), 0);
        List<direction> targetRefLevel(leaves.size(), direction(0));
        bool useNLayers(false);

        //- select boxes which need to be refined
        # ifdef USE_OMP
        # pragma omp parallel for if( leaves.size() > 1000 ) \
        reduction( + : nMarked) schedule(dynamic, 20)
        # endif
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if( oc.cubeType() & meshOctreeCubeBasic::OUTSIDE )
                continue;

            boundBox bb;
            oc.cubeBox(rootBox, bb.min(), bb.max());

            bool refine(false);
            forAll(refObjects, oI)
            {
                if( refObjects[oI].intersectsObject(bb) )
                {
                    # ifdef DEBUGSearch
                    Info << "Marking leaf " << leafI
                        << " for refinement" << endl;
                    # endif

                    if( oc.level() < refLevels[oI] )
                        refine = true;

                    if( refThickness[oI] > VSMALL )
                    {
                        const scalar cs = bb.max().x() - bb.min().x();
                        nLayers[leafI] =
                            Foam::max
                            (
                                nLayers[leafI],
                                Foam::max(label(refThickness[oI]/cs), 1)
                            );

                        targetRefLevel[leafI] =
                            Foam::max(targetRefLevel[leafI], refLevels[oI]);
                        useNLayers = true;
                    }
                }
            }

            if( refine )
            {
                refineCubes[leafI] = 1;
                ++nMarked;
            }
        }

        //- mark additional boxes for refinement to achieve
        //- correct refinement distance
        reduce(useNLayers, maxOp<label>());
        if( useNLayers )
        {
            nMarked +=
                octreeModifier.markAdditionalLayers
                (
                    refineCubes,
                    nLayers,
                    targetRefLevel
                );
        }

        //- refine boxes
        octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);

        # ifdef OCTREETiming
        const scalar refTime = omp_get_wtime();
        Info << "Time for refinement " << (refTime-startIter) << endl;
        # endif

        if( octree_.neiProcs().size() != 0 )
        {
            reduce(nMarked, sumOp<label>());
            if( nMarked )
            {
                octreeModifier.distributeLeavesToProcessors();

                # ifdef OCTREETiming
                const scalar distTime = omp_get_wtime();
                Info << "Time for distributing data to processors "
                << (distTime-refTime) << endl;
                # endif

                loadDistribution(false);

                # ifdef OCTREETiming
                Info << "Time for load distribution "
                << (omp_get_wtime()-distTime) << endl;
                # endif
            }
        }

    } while( nMarked != 0 );

    Info << "Finished refinement of boxes inside objects" << endl;

    //- set up inside-outside information
    createInsideOutsideInformation();
}

void meshOctreeCreator::refineBoxesIntersectingSurfaces()
{
    if( !meshDictPtr_ || !meshDictPtr_->found("surfaceMeshRefinement") )
    {
        return;
    }

    Info << "Refining boxes intersecting surface meshes" << endl;

    label nMarked;

    //- read surface meshes and calculate the refinement level for each
    //- surface mesh
    const dictionary& surfDict = meshDictPtr_->subDict("surfaceMeshRefinement");
    const wordList surfaces = surfDict.toc();
    PtrList<triSurf> surfaceMeshesPtr(surfaces.size());
    List<direction> refLevels(surfaces.size(), globalRefLevel_);
    scalarList refThickness(surfaces.size());

    //- load surface meshes into memory
    forAll(surfaceMeshesPtr, surfI)
    {
        const dictionary& dict = surfDict.subDict(surfaces[surfI]);

        const fileName fName(dict.lookup("surfaceFile"));

        surfaceMeshesPtr.set
        (
            surfI,
            new triSurf(fName)
        );

        direction addLevel(0);
        if( dict.found("cellSize") )
        {
            scalar s(readScalar(meshDictPtr_->lookup("maxCellSize")));

            const scalar cs = readScalar(dict.lookup("cellSize"));

            do
            {
                nMarked = 0;
                if( cs <= s * (1.+SMALL) )
                {
                    ++nMarked;
                    ++addLevel;
                }

                s /= 2.0;

            } while( nMarked != 0 );
        }
        else if( dict.found("additionalRefinementLevels") )
        {
            addLevel =
                readLabel(dict.lookup("additionalRefinementLevels"));
        }

        if( dict.found("refinementThickness") )
        {
            refThickness[surfI] =
                readScalar(dict.lookup("refinementThickness"));
        }

        //- set the refinement level for the current surface
        refLevels[surfI] += addLevel;
    }

    if( octree_.neiProcs().size() )
        forAll(refLevels, oI)
        {
            label l = refLevels[oI];
            reduce(l, maxOp<label>());
            refLevels[oI] = l;
        }

    //- start refining boxes intersecting triangles in each refinement surface
    const boundBox& rootBox = octree_.rootBox();
    const vector tol = SMALL * rootBox.span();
    meshOctreeModifier octreeModifier(octree_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();
    DynList<label> leavesInBox, intersectedLeaves;

    do
    {
        # ifdef OCTREETiming
        const scalar startIter = omp_get_wtime();
        # endif

        nMarked = 0;

        List<direction> refineCubes(leaves.size(), direction(0));
        labelList nLayers(leaves.size(), 0);
        List<direction> targetRefLevel(leaves.size(), direction(0));
        bool useNLayers(false);

        //- select boxes which need to be refined
        forAll(surfaceMeshesPtr, surfI)
        {
            const triSurf& surf = surfaceMeshesPtr[surfI];
            const pointField& points = surf.points();

            # ifdef USE_OMP
            # pragma omp parallel for \
            reduction( + : nMarked) schedule(dynamic, 10) \
            private(leavesInBox,intersectedLeaves)
            # endif
            forAll(surf, triI)
            {
                //- find the bounding box of the current triangle
                const labelledTri& tri = surf[triI];
                boundBox triBB(points[tri[0]], points[tri[0]]);
                for(label pI=1;pI<3;++pI)
                {
                    triBB.min() = Foam::min(triBB.min(), points[tri[pI]]);
                    triBB.max() = Foam::max(triBB.max(), points[tri[pI]]);
                }

                triBB.min() -= tol;
                triBB.max() += tol;

                //- find octree leaves inside the bounding box
                leavesInBox.clear();
                octree_.findLeavesContainedInBox(triBB, leavesInBox);

                //- check which of the leaves are intersected by the triangle
                intersectedLeaves.clear();
                forAll(leavesInBox, i)
                {
                    const label leafI = leavesInBox[i];

                    const meshOctreeCube& oc = *leaves[leafI];

                    if( oc.intersectsTriangleExact(surf, rootBox, triI) )
                    {
                        intersectedLeaves.append(leafI);

                        if( oc.level() < refLevels[surfI] )
                        {
                            # ifdef DEBUGSearch
                            Info << "Marking leaf " << leafI
                                << " with coordinates " << oc
                                << " for refinement" << endl;
                            # endif

                            if( !refineCubes[leafI] )
                                ++nMarked;
                            refineCubes[leafI] = 1;
                        }
                    }
                }

                if( refThickness[surfI] > VSMALL )
                {
                    useNLayers = true;

                    forAll(intersectedLeaves, i)
                    {
                        const label leafI = intersectedLeaves[i];
                        const meshOctreeCube& oc = *leaves[leafI];
                        const scalar cs = oc.size(rootBox);

                        nLayers[leafI] =
                            Foam::max
                            (
                                nLayers[leafI],
                                max(label(refThickness[surfI]/cs), 1)
                            );

                        targetRefLevel[leafI] =
                            Foam::max
                            (
                                targetRefLevel[leafI],
                                refLevels[surfI]
                            );
                    }
                }
            }
        }

        //- mark additional boxes for refinement to achieve
        //- correct refinement distance
        reduce(useNLayers, maxOp<label>());
        if( useNLayers )
        {
            nMarked +=
                octreeModifier.markAdditionalLayers
                (
                    refineCubes,
                    nLayers,
                    targetRefLevel
                );
        }

        //- refine boxes
        octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);

        # ifdef OCTREETiming
        const scalar refTime = omp_get_wtime();
        Info << "Time for refinement " << (refTime-startIter) << endl;
        # endif

        if( octree_.neiProcs().size() != 0 )
        {
            reduce(nMarked, sumOp<label>());
            if( nMarked )
            {
                octreeModifier.distributeLeavesToProcessors();

                # ifdef OCTREETiming
                const scalar distTime = omp_get_wtime();
                Info << "Time for distributing data to processors "
                << (distTime-refTime) << endl;
                # endif

                loadDistribution(false);

                # ifdef OCTREETiming
                Info << "Time for load distribution "
                << (omp_get_wtime()-distTime) << endl;
                # endif
            }
        }
    } while( nMarked != 0 );

    Info << "Finished refinement of boxes intersecting surface meshes" << endl;
}

void meshOctreeCreator::refineBoxesNearDataBoxes(const direction nLayers)
{
    # ifdef OCTREETiming
    const scalar startTime = omp_get_wtime();
    # endif

    const FixedList<meshOctreeCubeCoordinates, 26>& rp =
        octree_.regularityPositions();

    Info << "Refining boxes near DATA boxes" << endl;

    //- ensure one to one matching with unknown boxes
    meshOctreeModifier octreeModifier(octree_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    # ifdef DEBUGSearch
    forAll(leaves, leafI)
        Info << "Leaf " << leafI << " is " << *leaves[leafI]
            << " type " << label(leaves[leafI]->cubeType()) << endl;
    # endif

    List<direction> refineBox(leaves.size(), direction(0));

    labelHashSet transferCoordinates;
    LongList<meshOctreeCubeCoordinates> checkCoordinates;

    # ifdef USE_OMP
    # pragma omp parallel for if( leaves.size() > 1000 ) \
    schedule(dynamic, 20)
    # endif
    forAll(leaves, leafI)
    {
        if( leaves[leafI]->hasContainedElements() )
        {
            const meshOctreeCube& oc = *leaves[leafI];

            # ifdef DEBUGSearch
            Info << "Refining boxes near box " << leafI
                << " with coordinates " << oc << endl;
            # endif

            for(label k=0;k<26;++k)
            {
                const label neiLabel =
                    octree_.findLeafLabelForPosition
                    (
                        oc.coordinates() + rp[k]
                    );

                if( neiLabel == meshOctreeCube::OTHERPROC )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    {
                        if( !transferCoordinates.found(leafI) )
                        {
                            transferCoordinates.insert(leafI);
                            checkCoordinates.append(oc.coordinates());
                        }
                    }

                    continue;
                }

                if( neiLabel == -1 )
                    continue;

                const meshOctreeCube* nei = leaves[neiLabel];
                if( nei->level() == oc.level() )
                    continue;
                if( nei->cubeType() & meshOctreeCubeBasic::OUTSIDE )
                    continue;

                # ifdef DEBUGSearch
                Info << "Adding neighbour " << *nei << " of type"
                    << label(nei->cubeType()) << endl;
                # endif

                refineBox[nei->cubeLabel()] = 1;
            }
        }
    }

    if( octree_.neiProcs().size() )
    {
        LongList<meshOctreeCubeCoordinates> receivedCoordinates;
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            checkCoordinates,
            receivedCoordinates
        );

        # ifdef USE_OMP
        # pragma omp parallel for if( receivedCoordinates.size() > 1000 ) \
        schedule(dynamic, 20)
        # endif
        forAll(receivedCoordinates, ccI)
        {
            const meshOctreeCubeCoordinates& cc = receivedCoordinates[ccI];
            for(label k=0;k<26;++k)
            {
                const label neiLabel =
                    octree_.findLeafLabelForPosition(cc + rp[k]);

                if( neiLabel < 0 )
                    continue;

                const meshOctreeCube* nei = leaves[neiLabel];
                if( nei->level() == cc.level() )
                    continue;
                if( nei->cubeType() & meshOctreeCubeBasic::OUTSIDE )
                    continue;

                # ifdef DEBUGSearch
                Info << "Adding neighbour " << *nei << " of type"
                    << label(nei->cubeType()) << endl;
                # endif

                refineBox[nei->cubeLabel()] = 1;
            }
        }
    }

    for(direction i=1;i<nLayers;i++)
    {
        if( Pstream::parRun() )
        {
            checkCoordinates.clear();
            transferCoordinates.clear();
        }

        # ifdef USE_OMP
        # pragma omp parallel for if( leaves.size() > 1000 ) \
        schedule(dynamic, 20)
        # endif
        forAll(leaves, leafI)
        {
            if( refineBox[leafI] == i )
            {
                const meshOctreeCube& oc = *leaves[leafI];

                # ifdef DEBUGSearch
                Info << "Refining boxes near box " << leafI
                    << " with coordinates " << oc << endl;
                # endif

                for(label k=0;k<26;++k)
                {
                    const label neiLabel =
                        octree_.findLeafLabelForPosition
                        (
                            oc.coordinates() + rp[k]
                        );

                    if( neiLabel == meshOctreeCube::OTHERPROC )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        {
                            if( !transferCoordinates.found(leafI) )
                            {
                                transferCoordinates.insert(leafI);
                                checkCoordinates.append(oc.coordinates());
                            }
                        }

                        continue;
                    }

                    if( neiLabel == -1 )
                        continue;

                    const meshOctreeCube* nei = leaves[neiLabel];
                    if( nei->level() == oc.level() )
                        continue;
                    if( nei->cubeType() & meshOctreeCubeBasic::OUTSIDE )
                        continue;

                    # ifdef DEBUGSearch
                    Info << "Layer " << label(i) << endl;
                    Info << "Adding neighbour " << *nei << " of type"
                        << label(nei->cubeType()) << endl;
                    # endif

                    refineBox[nei->cubeLabel()] = i + 1;
                }
            }
        }

        if( octree_.neiProcs().size() )
        {
            LongList<meshOctreeCubeCoordinates> receivedCoordinates;
            octree_.exchangeRequestsWithNeighbourProcessors
            (
                checkCoordinates,
                receivedCoordinates
            );

            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoordinates.size() > 1000 ) \
            schedule(dynamic, 20)
            # endif
            forAll(receivedCoordinates, ccI)
            {
                const meshOctreeCubeCoordinates& cc = receivedCoordinates[ccI];

                for(label k=0;k<26;++k)
                {
                    const label neiLabel =
                        octree_.findLeafLabelForPosition(cc + rp[k]);

                    if( neiLabel < 0 )
                        continue;

                    const meshOctreeCube* nei = leaves[neiLabel];
                    if( nei->level() == cc.level() )
                        continue;
                    if( nei->cubeType() & meshOctreeCubeBasic::OUTSIDE )
                        continue;

                    # ifdef DEBUGSearch
                    Info << "Adding neighbour " << *nei << " of type"
                        << label(nei->cubeType()) << endl;
                    # endif

                    refineBox[nei->cubeLabel()] = i + 1;
                }
            }
        }
    }

    //- refine cubes
    octreeModifier.refineSelectedBoxes(refineBox, hexRefinement_);

    # ifdef OCTREETiming
    const scalar refTime = omp_get_wtime();
    Info << "Time for refinement " << (refTime-startTime) << endl;
    # endif

    if( Pstream::parRun() )
    {
        octreeModifier.distributeLeavesToProcessors();
        loadDistribution();
    }

    # ifdef OCTREETiming
    const scalar distTime = omp_get_wtime();
    Info << "Time for load balancing " << (distTime-refTime) << endl;
    # endif

    //- set up inside-outside information
    createInsideOutsideInformation();

    # ifdef OCTREETiming
    Info << "Time for calculation of inside/outside information "
        << (omp_get_wtime()-distTime) << endl;
    # endif

    Info << "Finished refining boxes near DATA boxes" << endl;
}

void meshOctreeCreator::refineBoxes
(
    const direction refLevel,
    const direction cubeType
)
{
    label nRefined;
    meshOctreeModifier octreeMod(octree_);

    do
    {
        # ifdef OCTREETiming
        const scalar startIter = omp_get_wtime();
        # endif

        nRefined = 0;

        const LongList<meshOctreeCube*>& leaves = octreeMod.leavesAccess();

        List<direction> refineCubes(leaves.size(), direction(0));

        # ifdef USE_OMP
        # pragma omp parallel for if( leaves.size() > 1000 ) \
        reduction(+ : nRefined) schedule(dynamic, 20)
        # endif
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if( (oc.cubeType() & cubeType) && (oc.level() < refLevel) )
            {
                refineCubes[leafI] = 1;
                ++nRefined;
            }
        }

        //- refine selected cubes
        octreeMod.refineSelectedBoxes(refineCubes, hexRefinement_);

        # ifdef OCTREETiming
        const scalar refTime = omp_get_wtime();
        Info << "Time for refinement " << (refTime-startIter) << endl;
        # endif

        if( Pstream::parRun() )
        {
            reduce(nRefined, sumOp<label>());
            if( nRefined )
            {
                octreeMod.distributeLeavesToProcessors();

                # ifdef OCTREETiming
                const scalar distTime = omp_get_wtime();
                Info << "Time for distributing data to processors "
                    << (distTime-refTime) << endl;
                # endif

                loadDistribution();

                # ifdef OCTREETiming
                Info << "Time for load distribution "
                    << (omp_get_wtime()-distTime) << endl;
                # endif
            }
        }

    } while( nRefined != 0 );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
