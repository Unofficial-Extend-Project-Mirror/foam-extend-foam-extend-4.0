/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

\*---------------------------------------------------------------------------*/

#include "removeFaces.H"
#include "polyMesh.H"
#include "directTopoChange.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "syncTools.H"
#include "OFstream.H"
#include "indirectPrimitivePatch.H"
#include "foamTime.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(removeFaces, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::removeFaces::changeCellRegion
(
    const label cellI,
    const label oldRegion,
    const label newRegion,
    labelList& cellRegion
) const
{
    if (cellRegion[cellI] == oldRegion)
    {
        cellRegion[cellI] = newRegion;

        // Go through all neighbouring cells
        const labelList& cCells = mesh_.cellCells()[cellI];

        forAll(cCells, i)
        {
            changeCellRegion(cCells[i], oldRegion, newRegion, cellRegion);
        }
    }
}


Foam::label Foam::removeFaces::changeFaceRegion
(
    const labelList& cellRegion,
    const boolList& removedFace,
    const labelList& nFacesPerEdge,
    const label faceI,
    const label newRegion,

    labelList& faceRegion
) const
{
    // Count number of changed faces
    label nChanged = 0;

    if (faceRegion[faceI] == -1 && !removedFace[faceI])
    {
        faceRegion[faceI] = newRegion;

        nChanged = 1;

        // Get mesh data
        const labelListList& meshFaceEdges = mesh_.faceEdges();
        const labelListList& meshEdgeFaces = mesh_.edgeFaces();

        // Step to neighbouring faces across edges that will get removed
        const labelList& fEdges = meshFaceEdges[faceI];

        forAll(fEdges, i)
        {
            const label& edgeI = fEdges[i];

            if (nFacesPerEdge[edgeI] >= 0 && nFacesPerEdge[edgeI] <= 2)
            {
                const labelList& eFaces = meshEdgeFaces[edgeI];

                forAll(eFaces, j)
                {
                    nChanged += changeFaceRegion
                    (
                        cellRegion,
                        removedFace,
                        nFacesPerEdge,
                        eFaces[j],
                        newRegion,

                        faceRegion
                    );
                }
            }
        }
    }

    return nChanged;
}


// Mark all faces affected in any way by
// - removal of cells
// - removal of faces
// - removal of edges
// - removal of points
Foam::Xfer<Foam::boolList> Foam::removeFaces::affectedFaces
(
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    const labelList& facesToRemove,
    const labelHashSet& edgesToRemove,
    const labelHashSet& pointsToRemove
) const
{
    // Create a marker field for affected faces
    boolList affectedFace(mesh_.nFaces(), false);

    // Get mesh cells
    const cellList& meshCells = mesh_.cells();

    // Mark faces affected by removal of cells
    forAll(cellRegion, cellI)
    {
        const label& region = cellRegion[cellI];

        // Bug fix: map coarse cell from all fine neighbours
        // All cells will be removed and replacew with a new one from
        // the central point.  HJ, 9/Sep/2019
        if (region != -1)
        {
            // Get this cell (list of cell faces) and mark all of its faces
            const labelList& cFaces = meshCells[cellI];

            forAll(cFaces, cFaceI)
            {
                affectedFace[cFaces[cFaceI]] = true;
            }
        }
    }

    // Mark faces affected by removal of face
    forAll(facesToRemove, i)
    {
         affectedFace[facesToRemove[i]] = true;
    }

    // Get edge faces
    const labelListList& meshEdgeFaces = mesh_.edgeFaces();

    // Mark faces affected by removal of edges
    forAllConstIter(labelHashSet, edgesToRemove, iter)
    {
        // Get all faces of this edge and mark them
        const labelList& eFaces = meshEdgeFaces[iter.key()];
        forAll(eFaces, eFaceI)
        {
            affectedFace[eFaces[eFaceI]] = true;
        }
    }

    // Get point faces
    const labelListList& meshPointFaces = mesh_.pointFaces();

    // Mark faces affected by removal of points
    forAllConstIter(labelHashSet, pointsToRemove, iter)
    {
        // Get all faces of this point and mark them
        const labelList& pFaces = meshPointFaces[iter.key()];
        forAll(pFaces, pFaceI)
        {
            affectedFace[pFaces[pFaceI]] = true;
        }
    }

    return xferMove(affectedFace);
}


void Foam::removeFaces::writeOBJ
(
    const indirectPrimitivePatch& fp,
    const fileName& fName
)
{
    OFstream str(fName);
    Pout<< "removeFaces::writeOBJ : Writing faces to file "
        << str.name() << endl;

    const pointField& localPoints = fp.localPoints();

    forAll(localPoints, i)
    {
        meshTools::writeOBJ(str, localPoints[i]);
    }

    const faceList& localFaces = fp.localFaces();

    forAll(localFaces, i)
    {
        const face& f = localFaces[i];

        str<< 'f';

        forAll(f, fp)
        {
            str<< token::SPACE << f[fp] + 1;
        }
        str<< nl;
    }
}


Foam::face Foam::removeFaces::filterFace
(
    const labelHashSet& pointsToRemove,
    const label faceI
) const
{
    // Get the face
    const face& f = mesh_.faces()[faceI];

    // Create a new face as plain list
    labelList newFace(f.size(), -1);

    label newFp = 0;

    forAll(f, fp)
    {
        // Get face vertex
        const label& vertI = f[fp];

        if (!pointsToRemove.found(vertI))
        {
            // Point not found, it won't be removed, append the vertex
            newFace[newFp++] = vertI;
        }
    }

    // Resize and return the face
    newFace.setSize(newFp);

    return face(newFace);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::removeFaces::removeFaces
(
    const polyMesh& mesh,
    const scalar minCos
)
:
    mesh_(mesh),
    minCos_(minCos)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::removeFaces::compatibleRemoves
(
    const labelList& facesToRemove,
    labelList& cellRegion,
    labelList& cellRegionMaster,
    labelList& newFacesToRemove
) const
{
    // Get necessary mesh data
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    const label nCells = mesh_.nCells();

    // Resize the lists and set elements to -1 by default
    cellRegion.setSize(nCells);
    cellRegion = -1;

    cellRegionMaster.setSize(nCells);
    cellRegionMaster = -1;

    // Count regions
    label nRegions = 0;

    // Loop through initial set of faces to remove
    forAll(facesToRemove, i)
    {
        // Get face index
        const label& faceI = facesToRemove[i];

        if (!mesh_.isInternalFace(faceI))
        {
            FatalErrorIn
            (
                "removeFaces::compatibleRemoves(const labelList&"
                ", labelList&, labelList&, labelList&)"
            )   << "Attempting to remove boundary face with face index: "
                << faceI << nl
                << "This is not allowed. Check facesToRemove list."
                << abort(FatalError);
        }

        // Get owner and neighbour data
        const label& own = faceOwner[faceI];
        const label& nei = faceNeighbour[faceI];

        // Get region data for owner and neighbour
        label region0 = cellRegion[own];
        label region1 = cellRegion[nei];

        if (region0 == -1)
        {
            // Region 0 (owner) is not set

            if (region1 == -1)
            {
                // Region 1 (neighbour) is also not set, create new region
                cellRegion[own] = nRegions;
                cellRegion[nei] = nRegions;

                // Make owner (lowest numbered!) the master of the region and
                // increment the number of regions
                cellRegionMaster[nRegions] = own;
                ++nRegions;
            }
            else
            {
                // Region 1 (neighbour) is set, add owner to neighbour region
                cellRegion[own] = region1;

                // See if owner becomes the master of the region (if its index
                // is lower than the current master of the region)
                cellRegionMaster[region1] = min(own, cellRegionMaster[region1]);
            }
        }
        else
        {
            // Region 0 (owner) is set

            if (region1 == -1)
            {
                // Neighbour region is not set, add neighbour to owner region
                cellRegion[nei] = region0;

                // Note: nei has higher index than own so neighbour can't be the
                // master of this region
            }
            else if (region0 != region1)
            {
                // Both have regions. Keep lowest numbered region and master
                label freedRegion = -1;
                label keptRegion = -1;

                if (region0 < region1)
                {
                    changeCellRegion
                    (
                        nei,
                        region1,    // old region
                        region0,    // new region
                        cellRegion
                    );

                    keptRegion = region0;
                    freedRegion = region1;
                }
                else if (region1 < region0)
                {
                    changeCellRegion
                    (
                        own,
                        region0,    // old region
                        region1,    // new region
                        cellRegion
                    );

                    keptRegion = region1;
                    freedRegion = region0;
                }

                label master0 = cellRegionMaster[region0];
                label master1 = cellRegionMaster[region1];

                cellRegionMaster[freedRegion] = -1;
                cellRegionMaster[keptRegion] = min(master0, master1);
            }
        }
    }

    // Set size
    cellRegionMaster.setSize(nRegions);

    // Note:
    // It is legal to have cellRegionMaster = -1 if the region has been created
    // and then abandoned because it has been merged with another region
    // HJ, 6/Sep/2019
    
    // Various checks, additional scope for clarity and memory management
    // - master is lowest numbered in any region
    // - regions have more than 1 cell
    {
        // Number of cells per region
        labelList nCells(cellRegionMaster.size(), 0);

        // Loop through cell regions
        forAll(cellRegion, cellI)
        {
            // Get region for this cell
            const label& r = cellRegion[cellI];

            if (r != -1)
            {
                // Region found, increment number of cells per region
                ++nCells[r];

                if (cellI < cellRegionMaster[r])
                {
                    FatalErrorInFunction
                        << "Not lowest numbered!  Cell: " << cellI
                        << " region: " << r
                        << " region master: " << cellRegionMaster[r]
                        << abort(FatalError);
                }
            }
        }

        // Loop through all regions
        forAll(nCells, regionI)
        {
            if (nCells[regionI] == 1)
            {
                FatalErrorInFunction
                    << "Region " << regionI
                    << " has only " << nCells[regionI] << " cell in it."
                    << abort(FatalError);
            }
        }
    }


    // Count number of used regions
    label nUsedRegions = 0;

    forAll(cellRegionMaster, i)
    {
        if (cellRegionMaster[i] != -1)
        {
            ++nUsedRegions;
        }
    }

    // Recreate facesToRemove to be consistent with the cellRegions
    dynamicLabelList allFacesToRemove(facesToRemove.size());

    // Get number of internal faces
    const label nInternalFaces = mesh_.nInternalFaces();

    // Loop through internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = faceOwner[faceI];
        const label& nei = faceNeighbour[faceI];

        if (cellRegion[own] != -1 && cellRegion[own] == cellRegion[nei])
        {
            // Both owner and neighbour of the face will become the same cell so
            // we can add this face to final list of faces to be removed
            allFacesToRemove.append(faceI);
        }
    }

    // Transfer dynamic list into the ordinary list
    newFacesToRemove.transfer(allFacesToRemove);

    return nUsedRegions;
}


// ************************************************************************* //
