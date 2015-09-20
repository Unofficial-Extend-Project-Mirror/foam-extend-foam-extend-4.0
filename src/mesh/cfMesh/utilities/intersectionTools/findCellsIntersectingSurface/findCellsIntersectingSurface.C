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

#include "findCellsIntersectingSurface.H"
#include "polyMeshGen.H"
#include "polyMeshGenAddressing.H"
#include "triSurf.H"
#include "boundBox.H"
#include "helperFunctions.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void findCellsIntersectingSurface::generateOctree(const triSurf& surf)
{
    octreePtr_ = new meshOctree(surf);

    meshOctreeCreator(*octreePtr_).createOctreeWithRefinedBoundary(15, 15);
}

void findCellsIntersectingSurface::findIntersectedCells()
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    const labelList& owner = mesh_.owner();

    const vectorField& faceCentres = mesh_.addressingData().faceCentres();
    const vectorField& cellCentres = mesh_.addressingData().cellCentres();

    meshOctreeModifier octreeModifier(*octreePtr_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    intersectedCells_.setSize(cells.size());
    facetsIntersectingCell_.setSize(cells.size());

    const triSurf& surf = octreePtr_->surface();
    const VRWGraph& pointFacets = surf.pointFacets();
    const pointField& sp = surf.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(cells, cellI)
    {
        bool intersected(false);

        const cell& c = cells[cellI];

        //- find the bounding box of the cell
        boundBox bb(cellCentres[cellI], cellCentres[cellI]);

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
            {
                bb.min() = Foam::min(bb.min(), points[f[pI]]);
                bb.max() = Foam::max(bb.max(), points[f[pI]]);
            }
        }

        const vector spanCorr = 0.01 * bb.span();
        bb.max() += spanCorr;
        bb.min() -= spanCorr;

        //- find surface triangles within the bounding box
        DynList<label> leavesInBox;
        octreePtr_->findLeavesContainedInBox(bb, leavesInBox);
        labelHashSet triangles;
        forAll(leavesInBox, boxI)
        {
            const meshOctreeCube& oc = *leaves[leavesInBox[boxI]];

            if( oc.hasContainedElements() )
            {
                const meshOctreeSlot& slot = *oc.slotPtr();
                const label ceI = oc.containedElements();

                forAllRow(slot.containedTriangles_, ceI, tI)
                    triangles.insert(slot.containedTriangles_(ceI, tI));
            }
        }

        //- remove triangles which do not intersect the bounding box
        labelHashSet reasonableCandidates;

        forAllConstIter(labelHashSet, triangles, tIter)
        {
            const labelledTri& tri = surf[tIter.key()];

            boundBox obb(sp[tri[0]], sp[tri[0]]);
            for(label i=1;i<3;++i)
            {
                const point& v = sp[tri[i]];
                obb.min() = Foam::min(obb.min(), v);
                obb.max() = Foam::max(obb.max(), v);
            }

            const vector spanTriCorr = SMALL * obb.span();
            obb.min() -= spanTriCorr;
            obb.max() += spanTriCorr;

            if( obb.overlaps(bb) )
                reasonableCandidates.insert(tIter.key());
        }

        triangles.transfer(reasonableCandidates);

        //- check if any of the surface vertices is contained within the cell
        labelHashSet nodes, facetsInCell;
        forAllConstIter(labelHashSet, triangles, tIter)
        {
            const labelledTri& tri = surf[tIter.key()];

            forAll(tri, i)
                nodes.insert(tri[i]);
        }

        //- check which surface nodes are within the cell
        forAllConstIter(labelHashSet, nodes, nIter)
        {
            const point& p = sp[nIter.key()];

            if( !bb.contains(p) )
                continue;

            bool foundIntersection(false);
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                if( owner[c[fI]] == cellI )
                {
                    forAll(f, pI)
                    {
                        tetrahedron<point, point> tet
                        (
                            points[f[pI]],
                            points[f.prevLabel(pI)],
                            faceCentres[c[fI]],
                            cellCentres[cellI]
                        );

                        if( help::pointInTetrahedron(p, tet) )
                        {
                            intersected = true;
                            forAllRow(pointFacets, nIter.key(), ptI)
                            {
                                facetsInCell.insert
                                (
                                    pointFacets(nIter.key(), ptI)
                                );
                            }

                            foundIntersection = true;
                            break;
                        }
                    }

                    if( foundIntersection )
                        break;
                }
                else
                {
                    forAll(f, pI)
                    {
                        tetrahedron<point, point> tet
                        (
                            points[f[pI]],
                            points[f.nextLabel(pI)],
                            faceCentres[c[fI]],
                            cellCentres[cellI]
                        );

                        if( help::pointInTetrahedron(p, tet) )
                        {
                            intersected = true;
                            forAllRow(pointFacets, nIter.key(), ptI)
                            {
                                facetsInCell.insert
                                (
                                    pointFacets(nIter.key(), ptI)
                                );
                            }

                            foundIntersection = true;
                            break;
                        }
                    }

                    if( foundIntersection )
                        break;
                }
            }
        }

        //- check if any triangle in the surface mesh
        //- intersects any of the cell's faces
        forAllConstIter(labelHashSet, triangles, tIter)
        {
            if( facetsInCell.found(tIter.key()) )
                continue;

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                const bool intersect =
                    help::doFaceAndTriangleIntersect
                    (
                        surf,
                        tIter.key(),
                        f,
                        points
                    );

                if( intersect )
                {
                    intersected = true;
                    facetsInCell.insert(tIter.key());
                    break;
                }
            }
        }

        //- store the results for this cell
        intersectedCells_[cellI] = intersected;
        if( intersected )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                facetsIntersectingCell_.setRow(cellI, facetsInCell.toc());
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

findCellsIntersectingSurface::findCellsIntersectingSurface
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    octreePtr_(const_cast<meshOctree*>(&octree)),
    octreeGenerated_(false),
    intersectedCells_(),
    facetsIntersectingCell_()
{
    findIntersectedCells();
}

findCellsIntersectingSurface::findCellsIntersectingSurface
(
    polyMeshGen& mesh,
    const triSurf& surface
)
:
    mesh_(mesh),
    octreePtr_(NULL),
    octreeGenerated_(true),
    intersectedCells_(),
    facetsIntersectingCell_()
{
    generateOctree(surface);

    findIntersectedCells();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

findCellsIntersectingSurface::~findCellsIntersectingSurface()
{
    if( octreeGenerated_ )
        deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const boolList& findCellsIntersectingSurface::intersectedCells() const
{
    return intersectedCells_;
}

const VRWGraph& findCellsIntersectingSurface::facetsIntersectingCells() const
{
    return facetsIntersectingCell_;
}

void findCellsIntersectingSurface::addIntersectedCellsToSubset
(
    const word subsetName
)
{
    const label setId = mesh_.addCellSubset(subsetName);

    forAll(intersectedCells_, cellI)
        if( intersectedCells_[cellI] )
            mesh_.addCellToSubset(setId, cellI);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
