/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "foamTime.H"
#include "meshOps.H"
#include "IOmanip.H"
#include "triFace.H"
#include "objectMap.H"
#include "faceSetAlgorithm.H"
#include "cellSetAlgorithm.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebugWithName(IOList<objectMap>, "objectMapIOList", 0);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Compute mapping weights for modified entities
void dynamicTopoFvMesh::computeMapping
(
    const scalar matchTol,
    const bool skipMapping,
    const bool mappingOutput,
    const label faceStart,
    const label faceSize,
    const label cellStart,
    const label cellSize
)
{
    // Convex-set algorithm for cells
    cellSetAlgorithm cellAlgorithm
    (
        (*this),
        oldPoints_,
        edges_,
        faces_,
        cells_,
        owner_,
        neighbour_
    );

    label nInconsistencies = 0;
    scalar maxFaceError = 0.0, maxCellError = 0.0;
    DynamicList<scalar> cellErrors(10), faceErrors(10);
    DynamicList<objectMap> failedCells(10), failedFaces(10);

    // Compute cell mapping
    for (label cellI = cellStart; cellI < (cellStart + cellSize); cellI++)
    {
        label cIndex = cellsFromCells_[cellI].index();
        labelList& masterObjects = cellsFromCells_[cellI].masterObjects();

        if (skipMapping)
        {
            // Dummy map from cell[0]
            masterObjects = labelList(1, 0);
            cellWeights_[cellI].setSize(1, 1.0);
            cellCentres_[cellI].setSize(1, vector::zero);
        }
        else
        {
            // Obtain weighting factors for this cell.
            cellAlgorithm.computeWeights
            (
                cIndex,
                0,
                cellParents_[cIndex],
                polyMesh::cellCells(),
                masterObjects,
                cellWeights_[cellI],
                cellCentres_[cellI]
            );

            // Add contributions from subMeshes, if any.
            computeCoupledWeights
            (
                cIndex,
                cellAlgorithm.dimension(),
                masterObjects,
                cellWeights_[cellI],
                cellCentres_[cellI]
            );

            // Compute error
            scalar error = mag(1.0 - sum(cellWeights_[cellI]));

            if (error > matchTol)
            {
                bool consistent = false;

                // Check whether any edges lie on boundary patches.
                // These cells can have relaxed weights to account
                // for mild convexity.
                const cell& cellToCheck = cells_[cIndex];

                if (is2D())
                {
                    const labelList& parents = cellParents_[cIndex];

                    forAll(parents, cI)
                    {
                        const cell& pCell = polyMesh::cells()[parents[cI]];

                        forAll(pCell, fI)
                        {
                            const face& pFace = polyMesh::faces()[pCell[fI]];

                            if (pFace.size() == 3)
                            {
                                continue;
                            }

                            label fP = boundaryMesh().whichPatch(pCell[fI]);

                            // Disregard processor patches
                            if (getNeighbourProcessor(fP) > -1)
                            {
                                continue;
                            }

                            if (fP > -1)
                            {
                                consistent = true;
                                break;
                            }
                        }

                        if (consistent)
                        {
                            break;
                        }
                    }
                }
                else
                {
                    forAll(cellToCheck, fI)
                    {
                        const labelList& fE = faceEdges_[cellToCheck[fI]];

                        forAll(fE, eI)
                        {
                            label eP = whichEdgePatch(fE[eI]);

                            // Disregard processor patches
                            if (getNeighbourProcessor(eP) > -1)
                            {
                                continue;
                            }

                            if (eP > -1)
                            {
                                consistent = true;
                                break;
                            }
                        }

                        if (consistent)
                        {
                            break;
                        }
                    }
                }

                if (!consistent)
                {
                    nInconsistencies++;

                    // Add to list
                    cellErrors.append(error);

                    // Accumulate error stats
                    maxCellError = Foam::max(maxCellError, error);

                    failedCells.append
                    (
                        objectMap
                        (
                            cIndex,
                            cellParents_[cIndex]
                        )
                    );
                }
            }
        }
    }

    // Convex-set algorithm for faces
    faceSetAlgorithm faceAlgorithm
    (
        (*this),
        oldPoints_,
        edges_,
        faces_,
        cells_,
        owner_,
        neighbour_
    );

    // Compute face mapping
    for (label faceI = faceStart; faceI < (faceStart + faceSize); faceI++)
    {
        label fIndex = facesFromFaces_[faceI].index();
        labelList& masterObjects = facesFromFaces_[faceI].masterObjects();

        label patchIndex = whichPatch(fIndex);
        label neiProc = getNeighbourProcessor(patchIndex);

        // Skip mapping for internal / processor faces.
        if (patchIndex == -1 || neiProc > -1)
        {
            // Set dummy masters, so that the conventional
            // faceMapper doesn't crash-and-burn
            masterObjects = labelList(1, 0);

            continue;
        }

        if (skipMapping)
        {
            // Dummy map from patch[0]
            masterObjects = labelList(1, 0);
            faceWeights_[faceI].setSize(1, 1.0);
            faceCentres_[faceI].setSize(1, vector::zero);
        }
        else
        {
            // Obtain weighting factors for this face.
            faceAlgorithm.computeWeights
            (
                fIndex,
                boundaryMesh()[patchIndex].start(),
                faceParents_[fIndex],
                boundaryMesh()[patchIndex].faceFaces(),
                masterObjects,
                faceWeights_[faceI],
                faceCentres_[faceI]
            );

            // Add contributions from subMeshes, if any.
            computeCoupledWeights
            (
                fIndex,
                faceAlgorithm.dimension(),
                masterObjects,
                faceWeights_[faceI],
                faceCentres_[faceI]
            );

            // Compute error
            scalar error = mag(1.0 - sum(faceWeights_[faceI]));

            if (error > matchTol)
            {
                bool consistent = false;

                // Check whether any edges lie on bounding curves.
                // These faces can have relaxed weights to account
                // for addressing into patches on the other side
                // of the curve.
                const labelList& fEdges = faceEdges_[fIndex];

                forAll(fEdges, eI)
                {
                    if (checkBoundingCurve(fEdges[eI]))
                    {
                        consistent = true;
                    }
                }

                if (!consistent)
                {
                    nInconsistencies++;

                    // Add to list
                    faceErrors.append(error);

                    // Accumulate error stats
                    maxFaceError = Foam::max(maxFaceError, error);

                    failedFaces.append
                    (
                        objectMap
                        (
                            fIndex,
                            faceParents_[fIndex]
                        )
                    );
                }
            }
        }
    }

    if (nInconsistencies)
    {
        Pout<< " Mapping errors: "
            << " max cell error: " << maxCellError
            << " max face error: " << maxFaceError
            << endl;

        if (debug || mappingOutput)
        {
            if (failedCells.size())
            {
                Pout<< " failedCells: " << endl;

                forAll(failedCells, cellI)
                {
                    label index = failedCells[cellI].index();

                    Pout<< "  Cell: " << index
                        << "  Error: " << cellErrors[cellI]
                        << endl;
                }
            }

            if (failedFaces.size())
            {
                Pout<< " failedFaces: " << endl;

                forAll(failedFaces, faceI)
                {
                    label index = failedFaces[faceI].index();
                    label patchIndex = whichPatch(index);

                    const polyBoundaryMesh& boundary = boundaryMesh();

                    Pout<< "  Face: " << index
                        << "  Patch: " << boundary[patchIndex].name()
                        << "  Error: " << faceErrors[faceI]
                        << endl;
                }
            }

            if (debug > 3 || mappingOutput)
            {
                // Prepare lists
                labelList objects;
                scalarField weights;
                vectorField centres;

                if (failedCells.size())
                {
                    forAll(failedCells, cellI)
                    {
                        label cIndex = failedCells[cellI].index();

                        cellAlgorithm.computeWeights
                        (
                            cIndex,
                            0,
                            cellParents_[cIndex],
                            polyMesh::cellCells(),
                            objects,
                            weights,
                            centres,
                            true
                        );

                        computeCoupledWeights
                        (
                            cIndex,
                            cellAlgorithm.dimension(),
                            objects,
                            weights,
                            centres,
                            true
                        );

                        // Clear lists
                        objects.clear();
                        weights.clear();
                        centres.clear();
                    }
                }

                if (failedFaces.size())
                {
                    forAll(failedFaces, faceI)
                    {
                        label fIndex = failedFaces[faceI].index();
                        label patchIndex = whichPatch(fIndex);

                        const polyBoundaryMesh& boundary = boundaryMesh();

                        faceAlgorithm.computeWeights
                        (
                            fIndex,
                            boundary[patchIndex].start(),
                            faceParents_[fIndex],
                            boundary[patchIndex].faceFaces(),
                            objects,
                            weights,
                            centres,
                            true
                        );

                        computeCoupledWeights
                        (
                            fIndex,
                            faceAlgorithm.dimension(),
                            objects,
                            weights,
                            centres,
                            true
                        );

                        // Clear lists
                        objects.clear();
                        weights.clear();
                        centres.clear();
                    }
                }
            }
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
    bool& mappingOutput = *(static_cast<bool*>(thread->operator()(2)));
    label& faceStart = *(static_cast<label*>(thread->operator()(3)));
    label& faceSize = *(static_cast<label*>(thread->operator()(4)));
    label& cellStart = *(static_cast<label*>(thread->operator()(5)));
    label& cellSize = *(static_cast<label*>(thread->operator()(6)));

    // Now calculate addressing
    mesh.computeMapping
    (
        matchTol,
        skipMapping,
        mappingOutput,
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
    bool skipMapping,
    bool mappingOutput
)
{
    label nThreads = threader_->getNumThreads();

    // If mapping is being skipped, issue a warning.
    if (skipMapping)
    {
        Info<< " *** Mapping is being skipped *** " << endl;
    }

    // Check if single-threaded
    if (nThreads == 1)
    {
        computeMapping
        (
            matchTol,
            skipMapping,
            mappingOutput,
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
        Pout<< " Mapping Faces: " << index[0] << nl
            << " Mapping Cells: " << index[1] << endl;
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
            Pout<< " Load starts: " << tStarts[indexI] << nl
                << " Load sizes: " << tSizes[indexI] << endl;
        }
    }

    // Set the argument list for each thread
    forAll(hdl, i)
    {
        // Size up the argument list
        hdl[i].setSize(7);

        // Set match tolerance
        hdl[i].set(0, &matchTol);

        // Set the skipMapping flag
        hdl[i].set(1, &skipMapping);

        // Set the mappingOutput flag
        hdl[i].set(2, &mappingOutput);

        // Set the start/size indices
        hdl[i].set(3, &(tStarts[0][i]));
        hdl[i].set(4, &(tSizes[0][i]));
        hdl[i].set(5, &(tStarts[1][i]));
        hdl[i].set(6, &(tSizes[1][i]));
    }

    // Prior to multi-threaded operation,
    // force calculation of demand-driven data.
    polyMesh::cells();
    primitiveMesh::cellCells();

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(boundary, patchI)
    {
        boundary[patchI].faceFaces();
    }

    // Force calculation of demand-driven data on subMeshes
    initCoupledWeights();

    // Execute threads in linear sequence
    executeThreads(identity(nThreads), hdl, &computeMappingThread);
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
            Pout<< "Inserting mapping cell: " << cIndex
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
    dynamicLabelList masterCells(5);

    forAll(mapCells, cellI)
    {
        if (mapCells[cellI] < 0)
        {
            continue;
        }

        if (mapCells[cellI] < nOldCells_)
        {
            if (findIndex(masterCells, mapCells[cellI]) == -1)
            {
                masterCells.append(mapCells[cellI]);
            }
        }

        if (cellParents_.found(mapCells[cellI]))
        {
            const labelList& nParents = cellParents_[mapCells[cellI]];

            forAll(nParents, cI)
            {
                if (findIndex(masterCells, nParents[cI]) == -1)
                {
                    masterCells.append(nParents[cI]);
                }
            }
        }
    }

    cellParents_.set(cIndex, masterCells);
}


// Set fill-in mapping information for a particular face
void dynamicTopoFvMesh::setFaceMapping
(
    const label fIndex,
    const labelList& mapFaces
)
{
    label patch = whichPatch(fIndex);
    label neiProc = getNeighbourProcessor(patch);

    if (debug > 3)
    {
        const polyBoundaryMesh& boundary = boundaryMesh();

        word pName;

        if (patch == -1)
        {
            pName = "Internal";
        }
        else
        if (patch < boundary.size())
        {
            pName = boundaryMesh()[patch].name();
        }
        else
        {
            pName = "New patch: " + Foam::name(patch);
        }

        Pout<< "Inserting mapping face: " << fIndex
            << " patch: " << pName
            << " mapFaces: " << mapFaces
            << " neiProc: "  << neiProc
            << endl;
    }

    bool foundError = false;

    // Check to ensure that internal faces are not mapped
    // from any faces in the mesh
    if (patch == -1 && mapFaces.size())
    {
        foundError = true;
    }

    // Check to ensure that mapFaces is non-empty for physical patches
    if (patch > -1 && mapFaces.empty() && neiProc == -1)
    {
        foundError = true;
    }

    if (foundError)
    {
        writeVTK("mapFace_" + Foam::name(fIndex), fIndex, 2);

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
            << " Patch: "
            << (patch > -1 ? boundaryMesh()[patch].name() : "Internal") << nl
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

    // For internal / processor faces, bail out
    if (patch == -1 || neiProc > -1)
    {
        return;
    }

    // Update face-parents information
    dynamicLabelList masterFaces(5);

    forAll(mapFaces, faceI)
    {
        if (mapFaces[faceI] < 0)
        {
            continue;
        }

        if (mapFaces[faceI] < nOldFaces_)
        {
            if (findIndex(masterFaces, mapFaces[faceI]) == -1)
            {
                masterFaces.append(mapFaces[faceI]);
            }
        }
        else
        if (faceParents_.found(mapFaces[faceI]))
        {
            const labelList& nParents = faceParents_[mapFaces[faceI]];

            forAll(nParents, fI)
            {
                if (findIndex(masterFaces, nParents[fI]) == -1)
                {
                    masterFaces.append(nParents[fI]);
                }
            }
        }
    }

    faceParents_.set(fIndex, masterFaces);
}


} // End namespace Foam

// ************************************************************************* //
