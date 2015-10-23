/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    private member of meshToMesh.
    Calculates mesh to mesh addressing pattern (for each cell from one mesh
    find the closest cell centre in the other mesh).

\*---------------------------------------------------------------------------*/

#include "meshToMesh.H"
#include "SubField.H"

#include "octree.H"
#include "octreeDataCell.H"
#include "octreeDataFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const Foam::debug::tolerancesSwitch
meshToMesh::cellCentreDistanceTol
(
    "meshToMeshCellCentreDistanceTol",
    1e-3
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void meshToMesh::calcAddressing()
{
    if (debug)
    {
        Info<< "meshToMesh::calculateAddressing() : "
            << "calculating mesh-to-mesh cell addressing" << endl;
    }

    // set reference to cells
    const cellList& fromCells = fromMesh_.cells();
    const pointField& fromPoints = fromMesh_.points();

    // In an attempt to preserve the efficiency of linear search and prevent
    // failure, a RESCUE mechanism will be set up. Here, we shall mark all
    // cells next to the solid boundaries. If such a cell is found as the
    // closest, the relationship between the origin and cell will be examined.
    // If the origin is outside the cell, a global n-squared search is
    // triggered.

    // SETTING UP RESCUE

    // visit all boundaries and mark the cell next to the boundary.

    if (debug)
    {
        Info<< "meshToMesh::calculateAddressing() : "
            << "Setting up rescue" << endl;
    }

    List<bool> boundaryCell(fromCells.size(), false);

    // set reference to boundary
    const polyPatchList& patchesFrom = fromMesh_.boundaryMesh();

    forAll (patchesFrom, patchI)
    {
        // get reference to cells next to the boundary
        const unallocLabelList& bCells = patchesFrom[patchI].faceCells();

        forAll (bCells, faceI)
        {
            boundaryCell[bCells[faceI]] = true;
        }
    }

    treeBoundBox meshBb(fromPoints);

    scalar typDim = meshBb.avgDim()/(2.0*cbrt(scalar(fromCells.size())));

    treeBoundBox shiftedBb
    (
        meshBb.min(),
        meshBb.max() + vector(typDim, typDim, typDim)
    );

    if (debug)
    {
        Info<< "\nMesh" << endl;
        Info<< "   bounding box           : " << meshBb << endl;
        Info<< "   bounding box (shifted) : " << shiftedBb << endl;
        Info<< "   typical dimension      :" << shiftedBb.typDim() << endl;
    }

    // Wrap indices and mesh information into helper object
    octreeDataCell shapes(fromMesh_);

    octree<octreeDataCell> oc
    (
        shiftedBb,  // overall bounding box
        shapes,     // all information needed to do checks on cells
        1,          // min. levels
        10.0,       // max. size of leaves
        10.0        // maximum ratio of cubes v.s. cells
    );

    if (debug)
    {
        oc.printStats(Info);
    }

    cellAddresses
    (
        cellAddressing_,
        toMesh_.cellCentres(),
        fromMesh_,
        boundaryCell,
        oc,
        true
    );

    forAll (toMesh_.boundaryMesh(), patchi)
    {
        const polyPatch& toPatch = toMesh_.boundaryMesh()[patchi];

        if (cuttingPatches_.found(toPatch.name()))
        {

            boundaryAddressing_[patchi].setSize(toPatch.size());

            cellAddresses
            (
                boundaryAddressing_[patchi],
                toPatch.faceCentres(),
                fromMesh_,
                boundaryCell,
                oc,
                false
            );
        }
        else if
        (
            patchMap_.found(toPatch.name())
         && fromMeshPatches_.found(patchMap_.find(toPatch.name())())
        )
        {
            const polyPatch& fromPatch = fromMesh_.boundaryMesh()
            [
                fromMeshPatches_.find(patchMap_.find(toPatch.name())())()
            ];

            if (fromPatch.empty())
            {
                WarningIn("meshToMesh::calcAddressing()")
                    << "Source patch " << fromPatch.name()
                    << " has no faces. Not performing mapping for it."
                    << endl;
                boundaryAddressing_[patchi] = -1;
            }
            else
            {
                treeBoundBox wallBb(fromPatch.localPoints());
                scalar typDim =
                    wallBb.avgDim()/(2.0*sqrt(scalar(fromPatch.size())));

                treeBoundBox shiftedBb
                (
                    wallBb.min(),
                    wallBb.max() + vector(typDim, typDim, typDim)
                );

                // Wrap data for octree into container
                octreeDataFace shapes(fromPatch);

                octree<octreeDataFace> oc
                (
                    shiftedBb,  // overall search domain
                    shapes,     // all information needed to do checks on cells
                    1,          // min levels
                    20.0,       // maximum ratio of cubes v.s. cells
                    2.0
                );


                const vectorField::subField centresToBoundary =
                    toPatch.faceCentres();

                boundaryAddressing_[patchi].setSize(toPatch.size());

                treeBoundBox tightest;
                scalar tightestDist;

                forAll(toPatch, toi)
                {
                    tightest = wallBb;                  // starting search bb
                    tightestDist = Foam::GREAT;        // starting max distance

                    boundaryAddressing_[patchi][toi] = oc.findNearest
                    (
                        centresToBoundary[toi],
                        tightest,
                        tightestDist
                    );
                }
            }
        }
    }

    if (debug)
    {
        Info<< "meshToMesh::calculateAddressing() : "
            << "finished calculating mesh-to-mesh acell ddressing" << endl;
    }
}


void meshToMesh::cellAddresses
(
    labelList& cellAddr,
    const pointField& points,
    const fvMesh& fromMesh,
    const List<bool>& boundaryCell,
    const octree<octreeDataCell>& oc,
    bool forceFind
) const
{
    label nCellsOutsideAddressing = 0;

    // The implemented search method is a simple neighbour array search.
    // It starts from a cell zero, searches its neighbours and finds one
    // which is nearer to the target point than the current position.
    // The location of the "current position" is reset to that cell and
    // search through the neighbours continues. The search is finished
    // when all the neighbours of the cell are farther from the target
    // point than the current cell

    // set curCell label to zero (start)
    register label curCell = 0;

    // set reference to cell to cell addressing
    const vectorField& centresFrom = fromMesh.cellCentres();
    const labelListList& cc = fromMesh.cellCells();

    forAll (points, toI)
    {
        scalar localTol = cellCentreDistanceTol();

        bool isBoundary = false;

        // pick up target position
        const vector& p = points[toI];

        // Set the sqr-distance
        scalar distSqr = magSqr(p - centresFrom[curCell]);

        bool closer;

        do
        {
            closer = false;

            // Set the current list of neighbouring cells
            const labelList& neighbours = cc[curCell];

            forAll (neighbours, nI)
            {
                scalar curDistSqr =
                    magSqr(p - centresFrom[neighbours[nI]]);

                // Search through all the neighbours.
                // If the cell is closer, reset current cell and distance
                if (curDistSqr < (1 - SMALL)*distSqr)
                {
                    curCell = neighbours[nI];
                    distSqr = curDistSqr;
                    closer = true;    // a closer neighbour has been found
                }
            }
        } while (closer);

        cellAddr[toI] = -1;

        // Check point is actually in the nearest cell
        if (fromMesh.pointInCell(p, curCell))
        {
            cellAddr[toI] = curCell;
        }
        else
        {
            // If curCell is a boundary cell then the point maybe either
            // outside the domain or in an other region of the doamin,
            //  either way use the octree search to find it.
            if (boundaryCell[curCell])
            {
                isBoundary = true;
                cellAddr[toI] = oc.find(p);
            }
            else
            {
                // If not on the boundary search the neighbours
                bool found = false;

                // set the current list of neighbouring cells
                const labelList& neighbours = cc[curCell];

                forAll (neighbours, nI)
                {
                    // search through all the neighbours.
                    // If point is in neighbour reset current cell
                    if (fromMesh.pointInCell(p, neighbours[nI]))
                    {
                        cellAddr[toI] = neighbours[nI];
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    // If still not found search the neighbour-neighbours

                    // set the current list of neighbouring cells
                    const labelList& neighbours = cc[curCell];

                    forAll (neighbours, nI)
                    {
                        // set the current list of neighbour-neighbouring cells
                        const labelList& nn = cc[neighbours[nI]];

                        forAll (nn, nI)
                        {
                            // search through all the neighbours.
                            // If point is in neighbour reset current cell
                            if (fromMesh.pointInCell(p, nn[nI]))
                            {
                                cellAddr[toI] = nn[nI];
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }

                if (!found)
                {
                    // Still not found so use the octree
                    cellAddr[toI] = oc.find(p);
                }
            }

            if (cellAddr[toI] < 0)
            {
                nCellsOutsideAddressing++;

                if (isBoundary && forceFind)
                {
                    // If still not found, get the closest cell within the
                    // specified tolerance

                    forAll(fromMesh.boundary(), patchi)
                    {
                        const fvPatch& patch = fromMesh.boundary()[patchi];

                        word name = patch.name();

                        label patchID =
                            toMesh_.boundaryMesh().findPatchID(name);

                        label sizePatch = 0;
                        if (patchID > -1)
                        {
                            sizePatch = toMesh_.boundary()[patchID].size();
                        }

                        if
                        (
                            sizePatch > 0
                        )
                        {
                            forAll(patch, facei)
                            {
                                label celli = patch.faceCells()[facei];

                                const vector& centre = fromMesh.C()[celli];
                                if (mag(points[toI] - centre) < localTol)
                                {
                                    localTol = mag(points[toI] - centre);
                                    cellAddr[toI] = celli;
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    if (nCellsOutsideAddressing > 0)
    {
        Info<< "Found " << nCellsOutsideAddressing
            << " cells outside of the addressing" << nl
            << "Cell addressing size = " << cellAddr.size() << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
