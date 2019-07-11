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

InClass
    domainDecomposition

Description
    Private member of domainDecomposition.
    Decomposes the mesh into bits

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "IOstreams.H"
#include "SLPtrList.H"
#include "boolList.H"
#include "cellList.H"
#include "primitiveMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::domainDecomposition::addInterProcessorBoundaryData
(
    const label& procI,
    const label& procJ,
    const label faceI,

    List<SLList<label> >& interProcBoundaries,
    List<SLList<SLList< label> > >& interProcBFaces
) const
{
    // Algorithm: Syncronously loop through the two lists containing all
    // neighbouring processor for a given processor and check whether this
    // processor is already in the list of neighbours. If the processor is
    // already there, simply append the face. If it's not there, append the
    // processor and the first face.

    // Get iterators to beginning of the two lists
    SLList<label>::iterator curInterProcBdrsIter =
        interProcBoundaries[procI].begin();

    SLList<SLList<label> >::iterator curInterProcBFacesIter =
        interProcBFaces[procI].begin();

    // Helper flag to distinguish when a particular neighbour processor has been
    // encountered
    bool interProcBouFound = false;

    // WARNING: Synchronous SLList iterators: Assuming interProcBoundaries and
    // interProcBFaces are always in sync as they should be
    for
    (
        ;

        curInterProcBdrsIter != interProcBoundaries[procI].end()
     && curInterProcBFacesIter != interProcBFaces[procI].end();

      ++curInterProcBdrsIter,
      ++curInterProcBFacesIter
    )
    {
        if (curInterProcBdrsIter() == procJ)
        {
            // Inter-processor boundary exists. Mark that we found this
            // neighbour processor
            interProcBouFound = true;

            // Append the face into the list
            curInterProcBFacesIter().append(faceI);

            // Break out of the loop
            break;
        }
    }

    if (!interProcBouFound)
    {
        // Inter-processor boundary does not exist and needs to be created

        // Append procJ to procI list
        interProcBoundaries[procI].append(procJ);

        // Append face to procI. Note: SLList construct taking a single
        // label creates a list with only that element
        interProcBFaces[procI].append(SLList<label>(faceI));
    }

    // Return whether the inter processor boundary has been found in the list
    return interProcBouFound;
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::domainDecomposition::decomposeMesh(const bool filterEmptyPatches)
{
    // Decide which cell goes to which processor
    distributeCells();

    // Distribute the cells according to the given processor label

    // Calculate the addressing information for the original mesh
    if (debug)
    {
        Pout<< "\nCalculating original mesh data" << endl;
    }

    // Set references to the original mesh
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Access all faces to grab the zones
    const faceList& allFaces = mesh_.allFaces();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // loop through the list of processor labels for the cell and add the
    // cell shape to the list of cells for the appropriate processor

    if (debug)
    {
        Pout<< "\nDistributing cells to processors" << endl;
    }

    // Memory management
    {
        List<SLList<label> > procCellList(nProcs_);

        forAll (cellToProc_, celli)
        {
            if (cellToProc_[celli] >= nProcs_)
            {
                FatalErrorIn("domainDecomposition::decomposeMesh()")
                    << "Impossible processor label " << cellToProc_[celli]
                    << "for cell " << celli
                    << abort(FatalError);
            }
            else
            {
                procCellList[cellToProc_[celli]].append(celli);
            }
        }

        // Convert linked lists into normal lists
        forAll (procCellList, procI)
        {
            procCellAddressing_[procI] = procCellList[procI];
        }
    }

    if (debug)
    {
        Pout << "\nDistributing faces to processors" << endl;
    }

    // Loop through internal faces and decide which processor they belong to
    // First visit all internal faces. If cells at both sides belong to the
    // same processor, the face is an internal face. If they are different,
    // it belongs to both processors.

    // Memory management
    {
        List<SLList<label> > procFaceList(nProcs_);

        forAll (neighbour, faceI)
        {
            if (cellToProc_[owner[faceI]] == cellToProc_[neighbour[faceI]])
            {
                // Face internal to processor
                procFaceList[cellToProc_[owner[faceI]]].append(faceI);
            }
        }

        // Record number of internal faces on each processor
        forAll (procFaceList, procI)
        {
            nInternalProcFaces_[procI] = procFaceList[procI].size();
        }

        // Detect inter-processor boundaries

        // Neighbour processor for each subdomain
        List<SLList<label> > interProcBoundaries(nProcs_);

        // Face labels belonging to each inter-processor boundary
        List<SLList<SLList<label> > > interProcBFaces(nProcs_);

        // Rewrite:
        // Handling of coupled patches is changed.  HJ, 11/Apr/2018

        // Prepare collection of patch faces

        // For all processors, set the size of start index and patch size
        // lists to the number of patches in the mesh
        forAll (procPatchSize_, procI)
        {
            procPatchSize_[procI].setSize(patches.size());
            procPatchStartIndex_[procI].setSize(patches.size());
        }

        forAll (patches, patchI)
        {
            // Reset size for all processors
            forAll (procPatchSize_, procI)
            {
                procPatchSize_[procI][patchI] = 0;
            }
        }


        // Algorithm:
        // When running the decomposition in parallel, it is assumed that
        // the result will be used in dynamic load balancing
        // For load balancing, the processor patches need to match in order for
        // the patch-to-patch matching to work properly on mesh reconstruction
        // Therefore, all processor patches need to be split into the matching
        // and non-matching part by examining the cellToProc_ data on
        // the neighbour side. This has been prepared in patchNbrCellToProc_
        // HJ, 11/Apr/2018

        // The correct order of dumping the faces is to dump the
        // slave processor patches first, then internal faces and then
        // master processor patches
        // HJ, 23/Apr/2018

        // Dump current processor patch faces into new processor patches
        // that will be created in decomposition when running in parallel
        forAll (patches, patchI)
        {
            // Check the processor patch for which neighbour data exists
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && !patchNbrCellToProc_[patchI].empty()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Note: DO ONLY SLAVE SIDE
                if (!procPatch.master())
                {
                    // Get patch start
                    const label patchStart = patches[patchI].start();

                    // Get faceCells
                    const labelList& fc = patches[patchI].faceCells();

                    // Get neighbour cellToProc addressing across the interface
                    const labelList& curNbrPtc = patchNbrCellToProc_[patchI];

                    forAll (fc, patchFaceI)
                    {
                        // Local owner proc is looked up using faceCells
                        const label ownerProc = cellToProc_[fc[patchFaceI]];

                        // Neighbour proc is looked up directly
                        const label neighbourProc = curNbrPtc[patchFaceI];

                        // Check change in processor type across the processor
                        // boundary

                        // If procI and procJ are the same, the processor face
                        // will be merged, meaning that it remains in the (old)
                        // processor patch. If procI and procJ are different,
                        // this will be a new processor boundary created from
                        // the existing processor face and added afterwards
                        if (ownerProc != neighbourProc)
                        {
                            // Insert inter-processor data for ownerProc
                            addInterProcessorBoundaryData
                            (
                                ownerProc,     // Processor to append to
                                neighbourProc, // Processor to append
                                patchStart + patchFaceI, // Face index to append

                                interProcBoundaries,
                                interProcBFaces
                            );
                        }

                        // Note: cannot insert regular faces here, because they
                        // are out of sequence.  HJ, 24/Apr/2018

                    } // End for all patch faces
                } // End if this is slave processor
            } // End if this is a processor patch
        } // End for all patches

        // Internal mesh faces
        forAll (neighbour, faceI)
        {
            const label& ownerProc = cellToProc_[owner[faceI]];
            const label& neighbourProc = cellToProc_[neighbour[faceI]];

            // Check whether we'll end up on the new processor boundary (if
            // ownerProc and neighbourProc are different)
            if (ownerProc != neighbourProc)
            {
                // Insert inter-processor data for ownerProc
                addInterProcessorBoundaryData
                (
                    ownerProc,     // Processor to append to
                    neighbourProc, // Processor to append
                    faceI,         // Face index to append

                    interProcBoundaries,
                    interProcBFaces
                );

                // Insert inter-processor data for neighbourProc
                // Note: ownerProc and neighbourProc swapped compared to the
                // call above
                addInterProcessorBoundaryData
                (
                    neighbourProc, // Processor to append to
                    ownerProc,     // Processor to append
                    faceI,         // Face index to append

                    interProcBoundaries,
                    interProcBFaces
                );
            }
        }

        // Dump current processor patch faces into new processor patches
        // that will be created in decomposition when running in parallel
        forAll (patches, patchI)
        {
            // Check the processor patch for which neighbour data exists
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && !patchNbrCellToProc_[patchI].empty()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Note: DO ONLY MASTER SIDE
                if (procPatch.master())
                {
                    // Get patch start
                    const label patchStart = patches[patchI].start();

                    // Get faceCells
                    const labelList& fc = patches[patchI].faceCells();

                    // Get neighbour cellToProc addressing across the interface
                    const labelList& curNbrPtc = patchNbrCellToProc_[patchI];

                    forAll (fc, patchFaceI)
                    {
                        // Local owner proc is looked up using faceCells
                        const label ownerProc = cellToProc_[fc[patchFaceI]];

                        // Neighbour proc is looked up directly
                        const label neighbourProc = curNbrPtc[patchFaceI];

                        // Check change in processor type across the processor
                        // boundary

                        // If ownerProc and neighbourProc are the same, the
                        // processor face will be merged, meaning that it
                        // remains in the (old) processor patch
                        // If ownerProc and neighbourProc are different,
                        // this will be a new processor boundary created from
                        // the existing processor face and added afterwards
                        if (ownerProc != neighbourProc)
                        {
                            // Insert inter-processor data for ownerProc
                            addInterProcessorBoundaryData
                            (
                                ownerProc,     // Processor to append to
                                neighbourProc, // Processor to append
                                patchStart + patchFaceI, // Face index to append

                                interProcBoundaries,
                                interProcBFaces
                            );
                        }

                    } // End for all patch faces
                } // End if this is slave processor
            } // End if this is a processor patch
        } // End for all patches

        // Loop through patches. For cyclic boundaries detect inter-processor
        // faces; for all other, add faces to the face list and remember start
        // and size of all patches.

        forAll (patches, patchI)
        {
            // New patch: record start index for all processors
            forAll (procPatchSize_, procI)
            {
                procPatchStartIndex_[procI][patchI] =
                    procFaceList[procI].size();
            }

            const label patchStart = patches[patchI].start();

            // Do normal patches. Note: processor patches have already been
            // partially done and need special treatment
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && !patchNbrCellToProc_[patchI].empty()
            )
            {
                // Only collect faces where the owner and neighbour processor
                // index are the same. If owner and neighbour processor index
                // are different, the face was already collected into a separate
                // patch. HJ, 23/Apr/2018.

                // Cast to processor patch
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Is this slave processor?
                const bool isSlave = !procPatch.master();

                // Get face cells on this side
                const unallocLabelList& fc = patches[patchI].faceCells();

                // Get face cells on the other side (communicated during
                // distributeCells() call)
                const unallocLabelList& nfc = patchNbrFaceCells_[patchI];

                // Get neighbour cellToProc addressing across the interface
                const labelList& curNbrPtc = patchNbrCellToProc_[patchI];

                // Note: in order to avoid upper triangular ordering errors for
                // neighbours after the reconstruction, we need to make sure
                // that the master processor inserts possibly multiple faces for
                // a single owner cell in the correct order. For each owner
                // cell with multiple faces, we need to collect neighbour cell
                // indices on the other side and possibly swap the order of
                // adding faces. The same "swapping" of insertion order needs to
                // happen on the slave side, but now we sort on the remote
                // (master) data and not on the local (slave) data as we do for
                // master processor. VV, 16/Feb/2019.

                // A map containing a list of patch faces and neighbour (owner)
                // cells for master (slave) processor, for each owner cell (for
                // master, neighbour cell for slave) (key). In the pair, first
                // entry is the list of patch faces and the second one is the
                // list of neigbhour (for master proc or owner for slave proc)
                // cells. We'll sort the list of neighbouring cells and use the
                // indices to add the patch faces in the correct order later
                // on. VV, 16/Feb/2019.
                Map<Pair<dynamicLabelList> > faceCellToNbrFaceCells
                (
                    // Reasonable size estimate (minimum = 10)
                    max(10, fc.size()/10)
                );

                forAll (fc, patchFaceI)
                {
                    // Get cell on this side: make a copy for swapping if this
                    // is slave processor later on
                    label ownCellI = fc[patchFaceI];

                    // Local owner proc is looked up using faceCells
                    const label ownerProc = cellToProc_[ownCellI];

                    // Neighbour proc is looked up directly
                    const label neighbourProc = curNbrPtc[patchFaceI];

                    // If the owner and neighbour processor index is the same,
                    // the face remains in the processor patch.
                    // In load balancing, it will be re-merged on reconstruction
                    // HJ, 23/Apr/2018.
                    if (ownerProc == neighbourProc)
                    {
                        // Get the face cell on the other side: make a copy for
                        // swapping if this is slave processor
                        label neiCellI = nfc[patchFaceI];

                        if (isSlave)
                        {
                            // This is slave, need to swap owner and neighbour
                            const label tmpCellI = ownCellI;
                            ownCellI = neiCellI;
                            neiCellI = tmpCellI;
                        }

                        // Add into the map
                        if (faceCellToNbrFaceCells.found(ownCellI))
                        {
                            // This owner cell has already been visited
                            Pair<dynamicLabelList>& tll =
                                faceCellToNbrFaceCells[ownCellI];

                            // Insert patch face index into the first list
                            tll.first().append(patchFaceI);

                            // Insert neighbour into the second list
                            tll.second().append(neiCellI);
                        }
                        else
                        {
                            // This owner cell has not been visited

                            // Create dynamicList containing the patch face
                            dynamicLabelList pfl(4);
                            pfl.append(patchFaceI);

                            // Create dynamicList containing the neighbour
                            dynamicLabelList nl(4);
                            nl.append(neiCellI);

                            // Insert the pair of lists into the map
                            faceCellToNbrFaceCells.insert
                            (
                                ownCellI,
                                Pair<dynamicLabelList>(pfl, nl)
                            );
                        }
                    } // End if owner and neighbour going to some processor
                } // End for all patch faces

                // Get sorted table of contents from the map to make sure that
                // we insert faces in the correct order (on both sides, since
                // for the slave processor, the map actually contains list of
                // patch faces and owner cells (on my, slave side) for all
                // neighbour cells (on the other, master side)
                const labelList sortedOwn = faceCellToNbrFaceCells.sortedToc();

                // Loop through ordered owner cells with faces on processor
                // patch
                forAll (sortedOwn, i)
                {
                    // Get current owner cell
                    const label& ownCellI = sortedOwn[i];

                    // The cell has not been handled yet and removed from
                    // the map. Get the pair of lists
                    Pair<dynamicLabelList>& patchFacesAndNeighbours =
                        faceCellToNbrFaceCells[ownCellI];

                    // Get the list of patch faces and list of neighbours
                    dynamicLabelList& pfl = patchFacesAndNeighbours.first();
                    dynamicLabelList& nl = patchFacesAndNeighbours.second();

                    // Create sorted neighbour list.
                    // Notes:
                    // 1. Transfer the contents of the dynamic list,
                    // 2. Sorted on construction
                    SortableList<label> sortedNbrs(nl.xfer());

                    // Get sorted indices for correct patch face insertion
                    // order
                    const labelList& sortedIndices = sortedNbrs.indices();

                    // Loop through sorted indices
                    forAll (sortedIndices, j)
                    {
                        // Get the sorted index into the patch face
                        const label& sI = sortedIndices[j];

                        // Get patch face index
                        const label& patchFaceI = pfl[sI];

                        // Get owner processor. Note: by definition/filtering
                        // above, owner and neighbour proc are the same. In
                        // order to avoid using cellToProc list for slave
                        // processor (where ownCellI is actually found on the
                        // other side), use curNbrPtc with patchFace index
                        const label& ownerProc = curNbrPtc[patchFaceI];

                        // Add the patch face into the list in the correct
                        // order for owner processor. Note: map contains
                        // only faces where owner proc = neighbour proc, so
                        // no need to double check
                        procFaceList[ownerProc].append(patchStart + patchFaceI);

                        // Increment the number of faces for this patch
                        procPatchSize_[ownerProc][patchI]++;
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                // Cyclic patch requires special treatment
                const polyPatch& cPatch = patches[patchI];

                const label cycOffset = cPatch.size()/2;

                // Set reference to faceCells for both patches
                const labelList::subList firstFaceCells
                (
                    cPatch.faceCells(),
                    cycOffset
                );

                const labelList::subList secondFaceCells
                (
                    cPatch.faceCells(),
                    cycOffset,
                    cycOffset
                );

                forAll (firstFaceCells, patchFaceI)
                {
                    // Get owner and neighbour processor labels
                    const label& ownerProc =
                        cellToProc_[firstFaceCells[patchFaceI]];
                    const label& neighbourProc =
                        cellToProc_[secondFaceCells[patchFaceI]];

                    // Check whether ownerProc and neighbourProc are different
                    if (ownerProc != neighbourProc)
                    {
                        // This face becomes an inter-processor boundary face
                        // inter - processor patch face found. Go through
                        // the list of inside boundaries for the owner
                        // processor and try to find this inter-processor
                        // patch.
                        cyclicParallel_ = true;

                        // Insert inter-processor data for ownerProc and return
                        // whether the neighbour was already present in the list
                        const bool ownerInterProcFound =
                            addInterProcessorBoundaryData
                            (
                                ownerProc,     // Processor to append to
                                neighbourProc, // Processor to append
                                
                                patchStart + patchFaceI, // Face index to append

                                interProcBoundaries,
                                interProcBFaces
                            );

                        // Insert inter-processor data for neighbourProc and
                        // return whether the neighbour was already present in
                        // the list
                        const bool neighbourInterProcFound =
                            addInterProcessorBoundaryData
                            (
                                neighbourProc, // Processor to append to
                                ownerProc,     // Processor to append
                                
                                // Face index with offset to append
                                patchStart + cycOffset + patchFaceI,

                                interProcBoundaries,
                                interProcBFaces
                            );

                        if (ownerInterProcFound && !neighbourInterProcFound)
                        {
                            FatalErrorIn
                            (
                                "domainDecomposition::decomposeMesh()"
                            )   << "Inconsistency in inter-processor "
                                << "boundary lists for processors "
                                << ownerProc << " and "
                                << neighbourProc
                                << " in cyclic boundary matching"
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // Add the face
                        procFaceList[ownerProc].append(patchStart + patchFaceI);

                        // Increment the number of faces for this patch
                        procPatchSize_[ownerProc][patchI]++;

                        // Note: I cannot add the other side of the cyclic
                        // boundary here because this would violate the order.
                        // They will be added in a separate loop below
                        // HJ, 15/Jan/2001
                    }
                }

                // Ordering in cyclic boundaries is important.
                // Add the other half of cyclic faces for cyclic boundaries
                // that remain on the processor
                forAll (secondFaceCells, patchFaceI)
                {
                    if
                    (
                        cellToProc_[firstFaceCells[patchFaceI]]
                     == cellToProc_[secondFaceCells[patchFaceI]]
                    )
                    {
                        // This cyclic face remains on the processor
                        const label ownerProc =
                            cellToProc_[firstFaceCells[patchFaceI]];

                        // Add the second face
                        procFaceList[ownerProc].append
                            (patchStart + cycOffset + patchFaceI);

                        // Increment the number of faces for this patch
                        procPatchSize_[ownerProc][patchI]++;
                    }
                }
            }
            else
            {
                // Normal patch. Add faces to processor where the cell
                // next to the face lives
                const unallocLabelList& fc = patches[patchI].faceCells();

                forAll (fc, patchFaceI)
                {
                    const label curProc = cellToProc_[fc[patchFaceI]];

                    // Add the face
                    procFaceList[curProc].append(patchStart + patchFaceI);

                    // Increment the number of faces for this patch
                    procPatchSize_[curProc][patchI]++;
                }
            }
        }

        // Face zone treatment.  HJ, 27/Mar/2009
        // Face zones identified as global will be present on all CPUs
        List<SLList<label> > procZoneFaceList(nProcs_);

        if (decompositionDict_.found("globalFaceZones"))
        {
            // Get faze zone names and face zones
            const wordList fzNames
            (
                decompositionDict_.lookup("globalFaceZones")
            );

            const faceZoneMesh& fz = mesh_.faceZones();

            forAll (fzNames, nameI)
            {
                const label zoneID = fz.findZoneID(fzNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn("domainDecomposition::decomposeMesh()")
                        << "Unknown global face zone " << fzNames[nameI]
                        << nl << "Valid face zones are" << fz.names()
                        << exit(FatalError);
                }

                if (debug)
                {
                    Pout<< "Preserving global face zone " << fzNames[nameI]
                        << endl;
                }

                const faceZone& curFz = fz[zoneID];

                // Go through all the faces in the zone.  If the owner of the
                // face equals to current processor, it has already been added;
                // otherwise, add the face to all processor face lists
                forAll (curFz, faceI)
                {
                    const label& curFaceID = curFz[faceI];

                    if (curFaceID < owner.size())
                    {
                        // This is an active mesh face, and it already belongs
                        // to one CPU.  Find out which and add it to the others
                        const label& curProc = cellToProc_[owner[curFaceID]];

                        forAll (procZoneFaceList, procI)
                        {
                            if (procI != curProc)
                            {
                                procZoneFaceList[procI].append(curFaceID);
                            }
                        }
                    }
                    else
                    {
                        // This is a stand-alone face, add it to all processors
                        forAll (procFaceList, procI)
                        {
                            procZoneFaceList[procI].append(curFaceID);
                        }
                    }
                }
            }
        }

        // Convert linked lists into normal lists
        // Add inter-processor boundaries and remember start indices
        forAll (procFaceList, procI)
        {
            // Get internal and regular boundary processor faces
            const SLList<label>& curProcFaces = procFaceList[procI];

            // Get reference to processor face addressing
            labelList& curProcFaceAddressing = procFaceAddressing_[procI];

            labelList& curProcNeighbourProcessors =
                procNeighbourProcessors_[procI];

            labelList& curProcProcessorPatchSize =
                procProcessorPatchSize_[procI];

            labelList& curProcProcessorPatchStartIndex =
                procProcessorPatchStartIndex_[procI];

            // Calculate the size
            label nFacesOnProcessor = curProcFaces.size();

            for
            (
                SLList<SLList<label> >::const_iterator curInterProcBFacesIter =
                    interProcBFaces[procI].cbegin();
                curInterProcBFacesIter != interProcBFaces[procI].cend();
              ++curInterProcBFacesIter
            )
            {
                nFacesOnProcessor += curInterProcBFacesIter().size();
            }

            // Add stand-alone global zone faces
            nFacesOnProcessor += procZoneFaceList[procI].size();

            curProcFaceAddressing.setSize(nFacesOnProcessor);

            // Fill in the list. Calculate turning index.
            // Turning index will be -1 only for some faces on processor
            // boundaries, i.e. the ones where the current processor ID
            // is in the cell which is a face neighbour.
            // Turning index is stored as the sign of the face addressing list

            label nFaces = 0;

            // Add internal and boundary faces
            // Remember to increment the index by one such that the
            // turning index works properly.  HJ, 5/Dec/2001
            for
            (
                SLList<label>::const_iterator curProcFacesIter =
                    curProcFaces.cbegin();
                curProcFacesIter != curProcFaces.cend();
              ++curProcFacesIter
            )
            {
                curProcFaceAddressing[nFaces] = curProcFacesIter() + 1;
                nFaces++;
            }

            // Add inter-processor boundary faces. At the beginning of each
            // patch, grab the patch start index and size

            curProcNeighbourProcessors.setSize
            (
                interProcBoundaries[procI].size()
            );

            curProcProcessorPatchSize.setSize
            (
                interProcBoundaries[procI].size()
            );

            curProcProcessorPatchStartIndex.setSize
            (
                interProcBoundaries[procI].size()
            );

            label nProcPatches = 0;

            SLList<label>::iterator curInterProcBdrsIter =
                interProcBoundaries[procI].begin();

            SLList<SLList<label> >::iterator curInterProcBFacesIter =
                interProcBFaces[procI].begin();

            for
            (
                ;

                curInterProcBdrsIter != interProcBoundaries[procI].end()
             && curInterProcBFacesIter != interProcBFaces[procI].end();

              ++curInterProcBdrsIter,
              ++curInterProcBFacesIter
            )
            {
                curProcNeighbourProcessors[nProcPatches] =
                    curInterProcBdrsIter();

                // Get start index for processor patch
                curProcProcessorPatchStartIndex[nProcPatches] = nFaces;

                label& curSize =
                    curProcProcessorPatchSize[nProcPatches];

                curSize = 0;

                // Add faces for this processor boundary

                for
                (
                    SLList<label>::iterator curFacesIter =
                        curInterProcBFacesIter().begin();
                    curFacesIter != curInterProcBFacesIter().end();
                  ++curFacesIter
                )
                {
                    // Add the face

                    // Remember to increment the index by one such that the
                    // turning index works properly.  HJ, 5/Dec/2001
                    if (cellToProc_[owner[curFacesIter()]] == procI)
                    {
                        curProcFaceAddressing[nFaces] = curFacesIter() + 1;
                    }
                    else
                    {
                        // Turning face
                        curProcFaceAddressing[nFaces] = -(curFacesIter() + 1);
                    }

                    // Increment the size
                    curSize++;

                    nFaces++;
                }

                nProcPatches++;
            }

            // Record number of live faces
            nLiveProcFaces_[procI] = nFaces;

            // Add stand-alone face zone faces
            const SLList<label>& curProcZoneFaces = procZoneFaceList[procI];

            for
            (
                SLList<label>::const_iterator curProcZoneFacesIter =
                    curProcZoneFaces.cbegin();
                curProcZoneFacesIter != curProcZoneFaces.cend();
              ++curProcZoneFacesIter
            )
            {
                curProcFaceAddressing[nFaces] = curProcZoneFacesIter() + 1;
                nFaces++;
            }
        } // End for all processors
    } // End of memory management

    if (debug)
    {
        Pout << "\nCalculating processor boundary addressing" << endl;
    }

    // For every patch of processor boundary, find the index of the original
    // patch. Mis-alignment is caused by the fact that patches with zero size
    // are omitted. For processor patches, set index to -1.
    // At the same time, filter the procPatchSize_ and procPatchStartIndex_
    // lists to exclude zero-size patches

    forAll (procPatchSize_, procI)
    {
        // Make a local copy of old lists
        const labelList oldPatchSizes = procPatchSize_[procI];

        const labelList oldPatchStarts = procPatchStartIndex_[procI];

        labelList& curPatchSizes = procPatchSize_[procI];

        labelList& curPatchStarts = procPatchStartIndex_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

        curBoundaryAddressing.setSize
        (
            oldPatchSizes.size()
          + curProcessorPatchSizes.size()
        );

        label nPatches = 0;

        forAll (oldPatchSizes, patchI)
        {
            // If filterEmptyPatches is set to true, or a patch is a
            // processor patch, remove it
            if
            (
                (filterEmptyPatches || isA<processorPolyPatch>(patches[patchI]))
             && oldPatchSizes[patchI] == 0
            )
            {
                // Patch filtered: do nothing
            }
            else
            {
                curBoundaryAddressing[nPatches] = patchI;

                curPatchSizes[nPatches] = oldPatchSizes[patchI];

                curPatchStarts[nPatches] = oldPatchStarts[patchI];

                nPatches++;
            }
        }

        // reset to the size of live patches
        curPatchSizes.setSize(nPatches);
        curPatchStarts.setSize(nPatches);

        forAll (curProcessorPatchSizes, procPatchI)
        {
            curBoundaryAddressing[nPatches] = -1;

            nPatches++;
        }

        curBoundaryAddressing.setSize(nPatches);
    }

    if (debug)
    {
        Pout << "\nDistributing points to processors" << endl;
    }

    // For every processor, loop through the list of faces for the processor.
    // For every face, loop through the list of points and mark the point as
    // used for the processor. Collect the list of used points for the
    // processor.

    // Record number of live points on each processor
    labelList nLivePoints(nProcs_, 0);

    forAll (procPointAddressing_, procI)
    {
        // Dimension list to all points in the mesh.  HJ, 27/Mar/2009
        boolList pointLabels(mesh_.allPoints().size(), false);

        // Get reference to list of used faces
        const labelList& procFaceLabels = procFaceAddressing_[procI];

        // Collect the used points
        labelList& procPointLabels = procPointAddressing_[procI];

        procPointLabels.setSize(pointLabels.size());

        // Two-pass algorithm:
        // First loop through live faces and record them in the points list
        // Second, visit all inactive zone faces and record the points

        label nUsedPoints = 0;

        // First pass: live faces

        for (label faceI = 0; faceI < nLiveProcFaces_[procI]; faceI++)
        {
            // Because of the turning index, some labels may be negative
            const labelList& facePoints =
                allFaces[mag(procFaceLabels[faceI]) - 1];

            forAll (facePoints, pointI)
            {
                // Mark the point as used
                pointLabels[facePoints[pointI]] = true;
            }
        }

        forAll (pointLabels, pointI)
        {
            if (pointLabels[pointI])
            {
                procPointLabels[nUsedPoints] = pointI;

                nUsedPoints++;
            }
        }

        // Record number of live points
        nLivePoints[procI] = nUsedPoints;

        // Second pass: zone faces

        // Reset point usage list
        boolList pointLabelsSecondPass(mesh_.allPoints().size(), false);

        for
        (
            label faceI = nLiveProcFaces_[procI];
            faceI < procFaceLabels.size();
            faceI++
        )
        {
            // Because of the turning index, some labels may be negative
            const labelList& facePoints =
                allFaces[mag(procFaceLabels[faceI]) - 1];

            forAll (facePoints, pointI)
            {
                // Mark the point as used
                if (!pointLabels[facePoints[pointI]])
                {
                    pointLabelsSecondPass[facePoints[pointI]] = true;
                }
            }
        }

        forAll (pointLabelsSecondPass, pointI)
        {
            if (pointLabelsSecondPass[pointI])
            {
                procPointLabels[nUsedPoints] = pointI;

                nUsedPoints++;
            }
        }

        // Reset the size of used points
        procPointLabels.setSize(nUsedPoints);
    }

    // Gather data about globally shared points

    // Memory management
    {
        // Dimension list to all points in the mesh.  HJ, 27/Mar/2009
        labelList pointsUsage(mesh_.allPoints().size(), 0);

        // Globally shared points are the ones used by more than 2 processors
        // Size the list approximately and gather the points
        labelHashSet gSharedPoints
        (
            min(100, mesh_.nPoints()/1000)
        );

        // Loop through all the processors and mark up points used by
        // processor boundaries.  When a point is used twice, it is a
        // globally shared point

        for (label procI = 0; procI < nProcs_; procI++)
        {
            // Get list of face labels
            const labelList& curFaceLabels = procFaceAddressing_[procI];

            // Get start of processor faces
            const labelList& curProcessorPatchStarts =
                procProcessorPatchStartIndex_[procI];

            const labelList& curProcessorPatchSizes =
                procProcessorPatchSize_[procI];

            // Reset the lookup list
            pointsUsage = 0;

            forAll (curProcessorPatchStarts, patchI)
            {
                const label curStart = curProcessorPatchStarts[patchI];
                const label curEnd = curStart + curProcessorPatchSizes[patchI];

                for
                (
                    label faceI = curStart;
                    faceI < curEnd;
                    faceI++
                )
                {
                    // Mark the original face as used
                    // Remember to decrement the index by one (turning index)
                    // HJ, 5/Dec/2001
                    const label curF = mag(curFaceLabels[faceI]) - 1;

                    const face& f = allFaces[curF];

                    forAll (f, pointI)
                    {
                        if (pointsUsage[f[pointI]] == 0)
                        {
                            // Point not previously used
                            pointsUsage[f[pointI]] = patchI + 1;
                        }
                        else if (pointsUsage[f[pointI]] != patchI + 1)
                        {
                            // Point used by some other patch = global point!
                            gSharedPoints.insert(f[pointI]);
                        }
                    }
                }
            }
        }

        // Grab the result from the hash list
        globallySharedPoints_ = gSharedPoints.toc();
        sort(globallySharedPoints_);
    }
}


// ************************************************************************* //
