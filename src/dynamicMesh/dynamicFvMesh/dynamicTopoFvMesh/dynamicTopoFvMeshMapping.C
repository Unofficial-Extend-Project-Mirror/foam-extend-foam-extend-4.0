/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    dynamicTopoFvMesh

Description
    Functions specific to conservative mapping

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "meshOps.H"
#include "IOmanip.H"
#include "triFace.H"
#include "objectMap.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Compute mapping weights for modified entities
void dynamicTopoFvMesh::computeMapping
(
    const scalar matchTol,
    const bool skipMapping,
    const label faceStart,
    const label faceSize,
    const label cellStart,
    const label cellSize
)
{
    // Compute cell mapping
    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        label cIndex = cellsFromCells_[cellI].index();

        if (skipMapping)
        {
            // Set empty mapping parameters
            const labelList& mo = cellParents_[cIndex];

            cellsFromCells_[cellI].masterObjects() = mo;
            cellWeights_[cellI].setSize(mo.size(), (1.0/(mo.size() + VSMALL)));
            cellCentres_[cellI].setSize(mo.size(), vector::zero);
        }
        else
        {
            // Compute master objects for inverse-distance weighting
            computeParents
            (
                cIndex,
                cellParents_[cIndex],
                polyMesh::cellCells(),
                polyMesh::cellCentres(),
                3,
                cellsFromCells_[cellI].masterObjects()
            );
        }
    }

    // Compute face mapping
    for (label faceI = faceStart; faceI < (faceStart + faceSize); faceI++)
    {
        label fIndex = facesFromFaces_[faceI].index();
        label patchIndex = whichPatch(fIndex);

        // Skip mapping for internal faces.
        if (patchIndex == -1)
        {
            // Set dummy masters, so that the conventional
            // faceMapper doesn't incur a seg-fault.
            facesFromFaces_[faceI].masterObjects() = labelList(1, 0);
            continue;
        }

        if (skipMapping)
        {
            // Set empty mapping parameters
            const labelList& mo = faceParents_[fIndex];

            facesFromFaces_[faceI].masterObjects() = mo;
            faceWeights_[faceI].setSize(mo.size(), (1.0/(mo.size() + VSMALL)));
            faceCentres_[faceI].setSize(mo.size(), vector::zero);
        }
        else
        {
            // Compute master objects for inverse-distance weighting
            computeParents
            (
                fIndex,
                faceParents_[fIndex],
                boundaryMesh()[patchIndex].faceFaces(),
                boundaryMesh()[patchIndex].faceCentres(),
                2,
                facesFromFaces_[faceI].masterObjects()
            );
        }
    }
}


// Static equivalent for multiThreading
void dynamicTopoFvMesh::computeMappingThread(void *argument)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Recast the pointers for the argument
    scalar& matchTol  = *(static_cast<scalar*>(thread->operator()(0)));
    bool& skipMapping = *(static_cast<bool*>(thread->operator()(1)));
    label& faceStart  = *(static_cast<label*>(thread->operator()(2)));
    label& faceSize   = *(static_cast<label*>(thread->operator()(3)));
    label& cellStart  = *(static_cast<label*>(thread->operator()(4)));
    label& cellSize   = *(static_cast<label*>(thread->operator()(5)));

    // Now calculate addressing
    mesh.computeMapping
    (
        matchTol,
        skipMapping,
        faceStart, faceSize,
        cellStart, cellSize
    );

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }
}


// Routine to invoke threaded mapping
void dynamicTopoFvMesh::threadedMapping
(
    scalar matchTol,
    bool skipMapping
)
{
    label nThreads = threader_->getNumThreads();

    // If mapping is being skipped, issue a warning.
    if (skipMapping)
    {
        Info << " *** Mapping is being skipped *** " << endl;
    }

    // Check if single-threaded
    if (nThreads == 1)
    {
        computeMapping
        (
            matchTol,
            skipMapping,
            0, facesFromFaces_.size(),
            0, cellsFromCells_.size()
        );

        return;
    }

    // Set one handler per thread
    PtrList<meshHandler> hdl(nThreads);

    forAll(hdl, i)
    {
        hdl.set(i, new meshHandler(*this, threader()));
    }

    // Simple load-balancing scheme
    FixedList<label, 2> index(-1);
    FixedList<labelList, 2> tStarts(labelList(nThreads, 0));
    FixedList<labelList, 2> tSizes(labelList(nThreads, 0));

    index[0] = facesFromFaces_.size();
    index[1] = cellsFromCells_.size();

    if (debug > 2)
    {
        Info << " Mapping Faces: " << index[0] << endl;
        Info << " Mapping Cells: " << index[1] << endl;
    }

    forAll(index, indexI)
    {
        label j = 0, total = 0;

        while (index[indexI]--)
        {
            tSizes[indexI][(j = tSizes[indexI].fcIndex(j))]++;
        }

        for (label i = 1; i < tStarts[indexI].size(); i++)
        {
            tStarts[indexI][i] = tSizes[indexI][i-1] + total;

            total += tSizes[indexI][i-1];
        }

        if (debug > 2)
        {
            Info << " Load starts: " << tStarts[indexI] << endl;
            Info << " Load sizes: " << tSizes[indexI] << endl;
        }
    }

    // Set the argument list for each thread
    forAll(hdl, i)
    {
        // Size up the argument list
        hdl[i].setSize(6);

        // Set match tolerance
        hdl[i].set(0, &matchTol);

        // Set the skipMapping flag
        hdl[i].set(1, &skipMapping);

        // Set the start/size indices
        hdl[i].set(2, &(tStarts[0][i]));
        hdl[i].set(3, &(tSizes[0][i]));
        hdl[i].set(4, &(tStarts[1][i]));
        hdl[i].set(5, &(tSizes[1][i]));
    }

    // Prior to multi-threaded operation,
    // force calculation of demand-driven data.
    polyMesh::cells();
    primitiveMesh::cellCells();
    primitiveMesh::cellCentres();

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(boundary, patchI)
    {
        boundary[patchI].faceFaces();
        boundary[patchI].faceCentres();
    }

    // Execute threads in linear sequence
    executeThreads(identity(nThreads), hdl, &computeMappingThread);
}


// Compute parents for inverse-distance weighting
void dynamicTopoFvMesh::computeParents
(
    const label index,
    const labelList& mapCandidates,
    const labelListList& oldNeighbourList,
    const vectorField& oldCentres,
    const label dimension,
    labelList& parents
) const
{
    if (parents.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "void dynamicTopoFvMesh::computeParents\n"
            "(\n"
            "    const label index,\n"
            "    const labelList& mapCandidates,\n"
            "    const labelListList& oldNeighbourList,\n"
            "    const vectorField& oldCentres,\n"
            "    const label dimension,\n"
            "    labelList& parents\n"
            ") const\n"
        )
            << " Addressing has already been calculated." << nl
            << " Index: " << index << nl
            << " Type: " << (dimension == 2 ? "Face" : "Cell") << nl
            << " mapCandidates: " << mapCandidates << nl
            << " Parents: " << parents << nl
            << abort(FatalError);
    }

    // Figure out the patch offset and centre
    label offset = -1;
    vector centre = vector::zero;

    if (dimension == 2)
    {
        offset = boundaryMesh()[whichPatch(index)].start();

        centre = faces_[index].centre(oldPoints_);
    }
    else
    if (dimension == 3)
    {
        offset = 0;
        scalar dummyVol = 0.0;

        meshOps::cellCentreAndVolume
        (
            index,
            oldPoints_,
            faces_,
            cells_,
            owner_,
            centre,
            dummyVol
        );
    }

    // Maintain a check-list
    labelHashSet checked;

    // Insert candidates first
    forAll(mapCandidates, indexI)
    {
        checked.insert(mapCandidates[indexI] - offset);
    }

    // Loop for three outward levels from candidates
    for (label i = 0; i < 3; i++)
    {
        // Fetch the set of candidates
        const labelList checkList = checked.toc();

        forAll(checkList, indexI)
        {
            const labelList& oldNei = oldNeighbourList[checkList[indexI]];

            forAll(oldNei, entityI)
            {
                if (!checked.found(oldNei[entityI]))
                {
                    checked.insert(oldNei[entityI]);
                }
            }
        }
    }

    // Loop through accumulated candidates and fetch the nearest one.
    scalar minDist = GREAT;
    label minLocation = -1;

    const labelList checkList = checked.toc();

    forAll(checkList, indexI)
    {
        scalar dist = magSqr(oldCentres[checkList[indexI]] - centre);

        if (dist < minDist)
        {
            minDist = dist;
            minLocation = checkList[indexI];
        }
    }

    // Set the final list
    const labelList& neiList = oldNeighbourList[minLocation];

    parents.setSize(neiList.size() + 1, -1);

    // Fill indices
    forAll(neiList, indexI)
    {
        parents[indexI] = neiList[indexI] + offset;
    }

    // Set final index
    parents[neiList.size()] = minLocation + offset;
}


// Set fill-in mapping information for a particular cell
void dynamicTopoFvMesh::setCellMapping
(
    const label cIndex,
    const labelList& mapCells,
    bool addEntry
)
{
    if (addEntry)
    {
        if (debug > 3)
        {
            Info << "Inserting mapping cell: " << cIndex << nl
                 << " mapCells: " << mapCells
                 << endl;
        }

        // Insert index into the list, and overwrite if necessary
        label index = -1;

        forAll(cellsFromCells_, indexI)
        {
            if (cellsFromCells_[indexI].index() == cIndex)
            {
                index = indexI;
                break;
            }
        }

        if (index == -1)
        {
            meshOps::sizeUpList
            (
                objectMap(cIndex, labelList(0)),
                cellsFromCells_
            );
        }
        else
        {
            cellsFromCells_[index].masterObjects() = labelList(0);
        }
    }

    // Update cell-parents information
    labelHashSet masterCells;

    forAll(mapCells, cellI)
    {
        if (mapCells[cellI] < 0)
        {
            continue;
        }

        if (mapCells[cellI] < nOldCells_)
        {
            masterCells.insert(mapCells[cellI]);
        }
        else
        if (cellParents_.found(mapCells[cellI]))
        {
            const labelList& nParents = cellParents_[mapCells[cellI]];

            forAll(nParents, cI)
            {
                masterCells.insert(nParents[cI]);
            }
        }
    }

    cellParents_.set(cIndex, masterCells.toc());
}


// Set fill-in mapping information for a particular face
void dynamicTopoFvMesh::setFaceMapping
(
    const label fIndex,
    const labelList& mapFaces
)
{
    label patch = whichPatch(fIndex);

    if (debug > 3)
    {
        Info << "Inserting mapping face: " << fIndex << nl
             << " patch: " << patch << nl
             << " mapFaces: " << mapFaces
             << endl;
    }

    bool foundError = false;

    // Check to ensure that internal faces are not mapped
    // from any faces in the mesh
    if (patch == -1 && mapFaces.size())
    {
        foundError = true;
    }

    // Check to ensure that boundary faces map
    // only from other faces on the same patch
    if (patch > -1 && mapFaces.empty())
    {
        foundError = true;
    }

    if (foundError)
    {
        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::setFaceMapping\n"
            "(\n"
            "    const label fIndex,\n"
            "    const labelList& mapFaces\n"
            ")"
        )
            << nl << " Incompatible mapping. " << nl
            << "  Possible reasons: " << nl
            << "    1. No mapping specified for a boundary face; " << nl
            << "    2. Mapping specified for an internal face, " << nl
            << "       when none was expected." << nl << nl
            << " Face: " << fIndex << nl
            << " Patch: " << patch << nl
            << " Owner: " << owner_[fIndex] << nl
            << " Neighbour: " << neighbour_[fIndex] << nl
            << " mapFaces: " << mapFaces << nl
            << abort(FatalError);
    }

    // Insert addressing into the list, and overwrite if necessary
    label index = -1;

    forAll(facesFromFaces_, indexI)
    {
        if (facesFromFaces_[indexI].index() == fIndex)
        {
            index = indexI;
            break;
        }
    }

    if (index == -1)
    {
        meshOps::sizeUpList
        (
            objectMap(fIndex, labelList(0)),
            facesFromFaces_
        );
    }
    else
    {
        facesFromFaces_[index].masterObjects() = labelList(0);
    }

    // For internal faces, set dummy maps / weights, and bail out
    if (patch == -1)
    {
        return;
    }

    // Update face-parents information
    labelHashSet masterFaces;

    forAll(mapFaces, faceI)
    {
        if (mapFaces[faceI] < 0)
        {
            continue;
        }

        if (mapFaces[faceI] < nOldFaces_)
        {
            masterFaces.insert(mapFaces[faceI]);
        }
        else
        if (faceParents_.found(mapFaces[faceI]))
        {
            const labelList& nParents = faceParents_[mapFaces[faceI]];

            forAll(nParents, fI)
            {
                masterFaces.insert(nParents[fI]);
            }
        }
    }

    faceParents_.set(fIndex, masterFaces.toc());
}


} // End namespace Foam

// ************************************************************************* //
