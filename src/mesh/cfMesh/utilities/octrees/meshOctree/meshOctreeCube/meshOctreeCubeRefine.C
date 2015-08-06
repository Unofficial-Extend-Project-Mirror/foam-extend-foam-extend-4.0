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

#include "meshOctreeCube.H"
#include "VRWGraph.H"
#include "triSurf.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCube::findContainedEdges
(
    const triSurf& surface,
    const boundBox& rootBox
)
{
    const VRWGraph& faceEdges = surface.facetEdges();
    const VRWGraph& edgeFaces = surface.edgeFacets();
    const edgeLongList& edges = surface.edges();
    const pointField& points = surface.points();

    const VRWGraph& containedElements = activeSlotPtr_->containedTriangles_;
    VRWGraph& containedEdges = activeSlotPtr_->containedEdges_;

    DynList<label> addedEdges;
    labelHashSet addEdge;
    forAllRow(containedElements, containedElementsLabel_, tI)
    {
        const label facetI = containedElements(containedElementsLabel_, tI);

        forAllRow(faceEdges, facetI, feI)
        {
            const label edgeI = faceEdges(facetI, feI);

            if( addEdge.found(edgeI) )
                continue;

            if( edgeFaces.sizeOfRow(edgeI) != 2 )
                continue;

            if(
                surface[edgeFaces(edgeI, 0)].region() !=
                surface[edgeFaces(edgeI, 1)].region()
            )
            {
                const edge& edg = edges[edgeI];
                const point& s = points[edg.start()];
                const point& e = points[edg.end()];
                if( intersectsLine(rootBox, s, e) )
                {
                    addEdge.insert(edgeI);
                    addedEdges.append(edgeI);
                }
            }
        }
    }

    if( addedEdges.size() != 0 )
    {
        containedEdgesLabel_ = containedEdges.size();
        containedEdges.appendList(addedEdges);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCube::refineCube2D
(
    const triSurf& surface,
    const boundBox& rootBox,
    meshOctreeSlot* slotPtr
)
{
    if( !slotPtr )
        slotPtr = activeSlotPtr_;

    # ifdef DEBUGSearch
    Info << "Refining the cube " << *this << endl;
    # endif

    //- set the cube label to -1
    cubeLabel_ = -1;

    //- create subCubes
    FixedList<meshOctreeCube*, 8> subCubes;
    forAll(subCubes, scI)
        subCubes[scI] = NULL;

    //- create new cubes in the Z-order fashion
    for(label scI=0;scI<4;++scI)
    {
        const label cubeI = slotPtr->cubes_.size();
        const meshOctreeCubeCoordinates cc = this->refineForPosition(scI);

        slotPtr->cubes_.append(cc);

        subCubes[scI] = &slotPtr->cubes_[cubeI];
        subCubes[scI]->activeSlotPtr_ = slotPtr;
        subCubes[scI]->setCubeType(this->cubeType());
        subCubes[scI]->setProcNo(this->procNo());
    }

    const label subCubesLabel = slotPtr->childCubes_.size();
    slotPtr->childCubes_.appendFixedList(subCubes);
    subCubesPtr_ = &slotPtr->childCubes_(subCubesLabel, 0);

    if( hasContainedElements() )
    {
        const VRWGraph& containedElements =
            activeSlotPtr_->containedTriangles_;

        # ifdef DEBUGSearch
        Info << "Distributing contained elements "
            << containedElements[containedElementsLabel_] << endl;
        # endif

        //- check if the subCube contain the element
        FixedList<DynList<label, 512>, 4> elementsInSubCubes;

        forAllRow(containedElements, containedElementsLabel_, tI)
        {
            const label elI = containedElements(containedElementsLabel_, tI);

            bool used(false);
            for(label scI=0;scI<4;++scI)
                if(
                    subCubes[scI]->intersectsTriangleExact
                    (
                        surface,
                        rootBox,
                        elI
                    )
                )
                {
                    used = true;
                    elementsInSubCubes[scI].append(elI);
                }

            if( !used )
            {
                Warning << "Triangle " << elI
                    << " is not transferred to the child cubes!" << endl;
            }
        }

        forAll(elementsInSubCubes, scI)
        {
            const DynList<label, 512>& elmts = elementsInSubCubes[scI];

            if( elmts.size() != 0 )
            {
                VRWGraph& ct = slotPtr->containedTriangles_;
                subCubes[scI]->containedElementsLabel_ = ct.size();
                ct.appendList(elmts);


                # ifdef DEBUGSearch
                Info << "Elements in leaf " << scI << " are "
                << ct[subCubes[scI]->containedElements()]
                    << endl;
                # endif
            }
        }

        //- find surface edges within the cube
        for(label scI=0;scI<4;++scI)
            if( subCubes[scI]->hasContainedElements() )
            {
                subCubes[scI]->findContainedEdges
                (
                    surface,
                    rootBox
                );
            }
            else if( subCubes[scI]->cubeType() & DATA )
            {
                subCubes[scI]->setCubeType(UNKNOWN);
            }
    }

    # ifdef DEBUGSearch
    for(label scI=0;scI<4;++scI)
        Info << "Refined cube " << scI << " is "
            << *subCubes[scI] << endl;
    # endif
}

void meshOctreeCube::refineCube
(
    const triSurf& surface,
    const boundBox& rootBox,
    meshOctreeSlot* slotPtr
)
{
    if( !slotPtr )
        slotPtr = activeSlotPtr_;

    # ifdef DEBUGSearch
    Info << "Refining the cube " << *this << endl;
    # endif

    //- set the cube label to -1
    cubeLabel_ = -1;

    //- create subCubes
    FixedList<meshOctreeCube*, 8> subCubes;

    //- create new cubes in the Z-order fashion
    for(label scI=0;scI<8;++scI)
    {
        const label cubeI = slotPtr->cubes_.size();
        slotPtr->cubes_.append(meshOctreeCube(this->refineForPosition(scI)));

        subCubes[scI] = &slotPtr->cubes_[cubeI];
        subCubes[scI]->activeSlotPtr_ = slotPtr;
        subCubes[scI]->setCubeType(this->cubeType());
        subCubes[scI]->setProcNo(this->procNo());
    }

    const label subCubesLabel = slotPtr->childCubes_.size();
    slotPtr->childCubes_.appendFixedList(subCubes);
    subCubesPtr_ = &slotPtr->childCubes_(subCubesLabel, 0);

    if( hasContainedElements() )
    {
        const VRWGraph& containedElements =
            activeSlotPtr_->containedTriangles_;

        # ifdef DEBUGSearch
        Info << "Distributing contained elements "
            << containedElements[containedElementsLabel_] << endl;
        # endif

        //- check if the subCube contain the element
        FixedList<DynList<label, 512>, 8> elementsInSubCubes;

        forAllRow(containedElements, containedElementsLabel_, tI)
        {
            const label elI = containedElements(containedElementsLabel_, tI);

            bool used(false);
            for(label scI=0;scI<8;++scI)
                if(
                    subCubes[scI]->intersectsTriangleExact
                    (
                        surface,
                        rootBox,
                        elI
                    )
                )
                {
                    used = true;
                    elementsInSubCubes[scI].append(elI);
                }

            if( !used )
            {
                Warning << "Triangle " << elI
                    << " is not transferred to the child cubes!" << endl;
            }
        }

        forAll(elementsInSubCubes, scI)
        {
            const DynList<label, 512>& elmts = elementsInSubCubes[scI];

            if( elmts.size() != 0 )
            {
                VRWGraph& ct = slotPtr->containedTriangles_;
                subCubes[scI]->containedElementsLabel_ = ct.size();
                ct.appendList(elmts);


                # ifdef DEBUGSearch
                Info << "Elements in leaf " << scI << " are "
                << ct[subCubes[scI]->containedElements()]
                    << endl;
                # endif
            }
        }

        //- find surface edges within the cube
        for(label scI=0;scI<8;++scI)
            if( subCubes[scI]->hasContainedElements() )
            {
                subCubes[scI]->findContainedEdges
                (
                    surface,
                    rootBox
                );
            }
            else if( subCubes[scI]->cubeType() & DATA )
            {
                subCubes[scI]->setCubeType(UNKNOWN);
            }
    }

    # ifdef DEBUGSearch
    for(label scI=0;scI<8;++scI)
        Info << "Refined cube " << scI << " is "
            << *subCubes[scI] << endl;
    # endif
}

void meshOctreeCube::refineMissingCube
(
    const triSurf& ts,
    const boundBox& rootBox,
    const label scI,
    meshOctreeSlot* slotPtr
)
{
    if( !slotPtr )
        slotPtr = activeSlotPtr_;

    if( !subCubesPtr_ )
    {
        FixedList<meshOctreeCube*, 8> sCubes;
        forAll(sCubes, i)
            sCubes[i] = NULL;

        const label subCubesLabel = slotPtr->childCubes_.size();

        slotPtr->childCubes_.appendFixedList(sCubes);
        subCubesPtr_ = &slotPtr->childCubes_(subCubesLabel, 0);
    }

    //- set the cube label to -1
    cubeLabel_ = -1;

    //- refine the cube for the selected position
    const label cubeI = slotPtr->cubes_.size();
    slotPtr->cubes_.append(meshOctreeCube(this->refineForPosition(scI)));
    subCubesPtr_[scI] = &slotPtr->cubes_[cubeI];

    subCubesPtr_[scI]->activeSlotPtr_ = slotPtr;
    subCubesPtr_[scI]->setCubeType(this->cubeType());
    subCubesPtr_[scI]->setProcNo(this->procNo());

    if( hasContainedElements() )
    {
        const VRWGraph& containedElements =
            activeSlotPtr_->containedTriangles_;
        DynList<label, 512> ce;

        forAllRow(containedElements, containedElementsLabel_, tI)
        {
            const label elI = containedElements(containedElementsLabel_, tI);
            if(
                subCubesPtr_[scI]->intersectsTriangleExact
                (
                    ts,
                    rootBox,
                    elI
                )
            )
            {
                ce.append(elI);
            }
        }

        if( ce.size() != 0 )
        {
            subCubesPtr_[scI]->containedElementsLabel_ =
                slotPtr->containedTriangles_.size();
            slotPtr->containedTriangles_.appendList(ce);

            subCubesPtr_[scI]->findContainedEdges
            (
                ts,
                rootBox
            );
        }
    }
}

void meshOctreeCube::refineMissingCube
(
    const label scI,
    const label elementsRowI,
    const label edgesRowI,
    meshOctreeSlot* slotPtr
)
{
    if( !slotPtr )
        slotPtr = activeSlotPtr_;

    if( subCubesPtr_ )
    {
        FixedList<meshOctreeCube*, 8> sCubes;
        forAll(sCubes, i)
            sCubes[i] = NULL;

        const label subCubesLabel = slotPtr->childCubes_.size();
        slotPtr->childCubes_.appendFixedList(sCubes);
        subCubesPtr_ = &slotPtr->childCubes_(subCubesLabel, 0);
    }

    //- set the cube label to -1
    cubeLabel_ = -1;

    //- refine the cube for the selected position
    const label cubeI = slotPtr->cubes_.size();
    slotPtr->cubes_.append(meshOctreeCube(this->refineForPosition(scI)));
    subCubesPtr_[scI] = &slotPtr->cubes_[cubeI];

    subCubesPtr_[scI]->activeSlotPtr_ = slotPtr;
    subCubesPtr_[scI]->setCubeType(this->cubeType());
    subCubesPtr_[scI]->setProcNo(this->procNo());

    //- set the contained elements and edges
    subCubesPtr_[scI]->containedElementsLabel_ = elementsRowI;
    subCubesPtr_[scI]->containedEdgesLabel_ = edgesRowI;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
