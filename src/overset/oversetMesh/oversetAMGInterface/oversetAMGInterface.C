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

#include "oversetAMGInterface.H"
#include "oversetMesh.H"
#include "IPstream.H"
#include "OPstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        AMGInterface,
        oversetAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetAMGInterface::initMap() const
{
    if (mapPtr_)
    {
        FatalErrorInFunction
            << "map already calculated"
            << abort(FatalError);
    }

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // Algorithm follows oversetMeshAddressing.C calcInterpolationMap
    // with changes.  Here, the information only exists on the acceptor
    // side.  HJ, 21/Sep/2019
    // Looping is performed over local acceptors

    const label nAcceptors = donorCellForAcceptor_.size();

    // Record number of donors from any processor.  This list will
    // be reduced
    labelListList nDonorsFromProcessorMap(Pstream::nProcs());

    // Collect donor cell index on the given processor
    labelListList donorFromProcessor(Pstream::nProcs());

    // Create construct map: what I am receiving from each processor
    labelListList constructMap(Pstream::nProcs());

    // Create list containing number of donors A processor is sending to other
    // processors.
    forAll (nDonorsFromProcessorMap, procI)
    {
        nDonorsFromProcessorMap[procI].setSize(Pstream::nProcs(), 0);

        donorFromProcessor[procI].setSize(nAcceptors);
        constructMap[procI].setSize(nAcceptors);
    }

    // Get the nDonorsFromProcessorMap for my processor
    labelList& numberOfLocalDonorsFromProcs =
        nDonorsFromProcessorMap[Pstream::myProcNo()];

    labelList nConstructMap(Pstream::nProcs(), 0);

    forAll (donorCellForAcceptor_, accI)
    {
        // Record donor cell from donor proc
        donorFromProcessor[donorProcForAcceptor_[accI]]
            [numberOfLocalDonorsFromProcs[donorProcForAcceptor_[accI]]] =
            donorProcIndexForAcceptor_[accI];

        // Record that the processor is giving data
        constructMap[donorProcForAcceptor_[accI]]
            [numberOfLocalDonorsFromProcs[donorProcForAcceptor_[accI]]] =
            accI;

        // Increment number of donors from this proc
        numberOfLocalDonorsFromProcs[donorProcForAcceptor_[accI]]++;
    }

    // Resize lists
    forAll (donorFromProcessor, procI)
    {
        donorFromProcessor[procI].setSize
        (
            numberOfLocalDonorsFromProcs[procI]
        );

        constructMap[procI].setSize
        (
            numberOfLocalDonorsFromProcs[procI]
        );
    }

    // Reduce the number of cells.  This will allow the receiving
    // processor to know it will receive the sendMap
    Pstream::gatherList(nDonorsFromProcessorMap);
    Pstream::scatterList(nDonorsFromProcessorMap);

    // Send the donor cell indices to other processor
    // Lookup MY ROW
    forAll (numberOfLocalDonorsFromProcs, nbrProcI)
    {
        // If there are entries, send them across
        if (numberOfLocalDonorsFromProcs[nbrProcI] > 0)
        {
            OPstream toNbrProc
            (
                Pstream::blocking,
                nbrProcI
            );

            toNbrProc << donorFromProcessor[nbrProcI];
        }
    }

    // Create sendMap.  All data will be received from nbr
    // processors
    labelListList sendMap(Pstream::nProcs());

    // Receive data from the other side.  Loop through second
    // index
    forAll (nDonorsFromProcessorMap, nbrProcI)
    {
        if (nDonorsFromProcessorMap[nbrProcI][Pstream::myProcNo()] > 0)
        {
            // There is data to receive
            IPstream fromNbrProc
            (
                Pstream::blocking,
                nbrProcI
            );

            // Insert into sendMap for nbrProc
            sendMap[nbrProcI] = labelList(fromNbrProc);
        }
    }

    mapPtr_ = new mapDistribute
    (
        acceptorCells().size(),
        sendMap,
        constructMap
    );

    // // Collect donor cells.  Use sendMap and unpack it
    label nDonorCells = 0;

    forAll (sendMap, procI)
    {
        nDonorCells += sendMap[procI].size();
    }

    // Set donor cells and donor proc
    donorCells_.setSize(nDonorCells);
    donorCellsProc_.setSize(nDonorCells);

    // Reset counter
    nDonorCells = 0;

    // Collect donor cells
    forAll (sendMap, procI)
    {
        const labelList& curSend = sendMap[procI];

        forAll (curSend, csI)
        {
            donorCells_[nDonorCells] = curSend[csI];
            donorCellsProc_[nDonorCells] = procI;

            nDonorCells++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetAMGInterface::oversetAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    AMGInterface(lduMesh),
    fineOversetInterface_(refCast<const oversetLduInterface>(fineInterface)),
    mapPtr_(nullptr)
{
    // Info<< "Constructing overset interface" << nl
    //     << "local restrict: " << localRestrictAddressing.size() << nl
    //     << "neighbour restrict: " << neighbourRestrictAddressing.size()
    //     << endl;

    // Size of coarse interface is equal to max cluster size plus one
    interfaceSize_ = max(localRestrictAddressing) + 1;

    // Note: cannot deal with weighted interpolation because access
    // to interpolation weights is coded under the name of the field being
    // solved for, which is not available for AMG coarsening
    // The coarsening across overset interface shall be done using only the
    // master acceptor and a single weight
    // HJ, 11/Sep/2019

    // Initialise fine map
    // Note: the mapDistribute contains a waitRequests call which cannot
    // be limited to a comm.  Therefore, all processors need to calculate
    // the mapDistribute schedule before escaping the constructor,
    // even if there are no overset cells available.
    // HJ, 18/Oct/2016
    if (Pstream::parRun())
    {
        fineOversetInterface_.map().schedule();
    }

    // Coarsen local donor cells and renumber.  Record local donor index
    // and communicate across

    // Note
    // For serial runs donors are calculated here
    // For parallel runs, they will be over-ridden in initMap
    // There are 2 ways of calculating donors
    // - go to fine donors and pick up restrict addressing
    // - take sendMap and unpack.  The coarse donor will have the same value
    //   but the donor that provides data to multiple neighbouring processors
    //   will be recorded multiple times.  However, its restrict index is
    //   the same.
    // Reconsider.  HJ, 23/Sep/2019
    
    // Get fine donors
    const labelList& fineDonors = fineOversetInterface_.donorCells();

    // Get fine donors processor (to which donor is sent)
    const labelList& fineDonorsProc = fineOversetInterface_.donorCellsProc();

    // Memory management
    {
        // Note: donor collection must be done per-processor: a single
        // donor may need to be passed to multiple acceptors on
        // multiple processors.  HJ, 23/Sep/2019

        // Into the hash table, using the coarse donor, record coarse processor
        // from fine donor.  They must be the same
        List<labelHashSet> coarseDonorProcMap(Pstream::nProcs());

        // For all fine donors, find and collect restrict addressing.
        // This will be the coarse donor
        forAll (fineDonors, fdI)
        {
            coarseDonorProcMap[fineDonorsProc[fdI]].insert
            (
                localRestrictAddressing[fineDonors[fdI]]
            );
        }

        // Collect coarse donor cells

        // Count donor cells
        label nDonors = 0;

        forAll (coarseDonorProcMap, procI)
        {
            nDonors+= coarseDonorProcMap[procI].size();
        }

        donorCells_.setSize(nDonors);
        donorCellsProc_.setSize(nDonors);

        // Reset counter
        nDonors = 0;

        forAll (coarseDonorProcMap, procI)
        {
            const labelList& curDonors = coarseDonorProcMap[procI].sortedToc();

            forAll (curDonors, cdI)
            {
                donorCells_[nDonors] = curDonors[cdI];
                donorCellsProc_[nDonors] = procI;

                nDonors++;
            }
        }
    }

    // Having established coarse donor list (in ascending order),
    // for each fine donor find the position of the coarse donor in the
    // list

    // I will know coarse donor index and I need its location in the list
    labelList coarseDonorIndex(max(donorCells_) + 1);

    forAll (donorCells_, dcI)
    {
        coarseDonorIndex[donorCells_[dcI]] = dcI;
    }

    // Make a list across fine level, where each fine donor records
    // coarse donor index
    labelList acceptorCellDonorIndex
    (
        fineOversetInterface_.interfaceSize(),
        -1
    );

    // Fill the list: for all fine donors, look up coarse donor,
    // look up coarse donor index and put it into acceptorCellDonorIndex
    forAll (fineDonors, fdI)
    {
        // Ger coarse donor index
        const label cdi = localRestrictAddressing[fineDonors[fdI]];

        // In the acceptor cell donor list, put coarse donor index
        acceptorCellDonorIndex[fineDonors[fdI]] = coarseDonorIndex[cdi];
    }

    // Communicate coarse donor
    fineOversetInterface_.initInternalFieldTransfer
    (
        Pstream::blocking,
        acceptorCellDonorIndex
    );

    acceptorCellDonorIndex =
        fineOversetInterface_.internalFieldTransfer
        (
            Pstream::blocking,
            acceptorCellDonorIndex
        );

    // Donor processor info needs to be communicated.
    // Use fine level comms
    labelList acceptorCellDonorProc
    (
        fineOversetInterface_.interfaceSize(),
        Pstream::myProcNo()
    );

    // Communicate donor proc
    fineOversetInterface_.initInternalFieldTransfer
    (
        Pstream::blocking,
        acceptorCellDonorProc
    );

    acceptorCellDonorProc =
        fineOversetInterface_.internalFieldTransfer
        (
            Pstream::blocking,
            acceptorCellDonorProc
        );

    // Note: Guessing size of HashTable to fine interface size

    // Coded neighbour index. Note: using long int to simplify encoding
    // HJ, 1/Aug/2016
    HashTable<DynamicList<long, 4>, long, Hash<long> > neighboursTable
    (
        Foam::max(128, fineOversetInterface_.interfaceSize()/4)
    );

    // Neighbour processor index
    HashTable<DynamicList<label, 4>, label, Hash<label> > nbrsProcTable
    (
        Foam::max(128, fineOversetInterface_.interfaceSize()/4)
    );

    // Neighbour processor donor index
    HashTable<DynamicList<label, 4>, label, Hash<label> > nbrsProcDonorIdTable
    (
        Foam::max(128, fineOversetInterface_.interfaceSize()/4)
    );

    // Neighbour face-faces addressing for a face with split neighbours
    HashTable<DynamicList<DynamicList<label, 4>, 4>, label, Hash<label> >
    faceFaceTable
    (
        Foam::max(128, fineOversetInterface_.interfaceSize()/4)
    );

    // Count the number of coarse faces
    label nCoarseFaces = 0;

    // Count the number of agglomeration pairs
    label nAgglomPairs = 0;

    // On the fine level, addressing is obtained from oversetMesh
    if (fineOversetInterface_.fineLevel())
    {
        // Loop over all fringe faces
        // Follow the pattern in oversetFvPatchField::initInterfaceMatrixUpdate

        // Get mesh owner-neighbour addressing to visit cells around fringe
        // faces
        const unallocLabelList& own =
            fineOversetInterface_.ldu().lduAddr().lowerAddr();

        const unallocLabelList& nei =
            fineOversetInterface_.ldu().lduAddr().upperAddr();

        const labelList& fringeFaces =
            fineOversetInterface_.overset().fringeFaces();

        const labelList& fringeFaceCells =
            fineOversetInterface_.overset().fringeFaceCells();

        const boolList& fringeFaceFlips =
            fineOversetInterface_.overset().fringeFaceFlips();

        label curMasterProc, curSlaveProc, curSlaveProcDonorId;
        long curMaster, curSlave;

        // On the fine level, loop through fringe faces:
        // Overset interfaces are between internal cells
        forAll (fringeFaces, ffI)
        {
            const label& curFace = fringeFaces[ffI];

            if (curFace < own.size())
            {
                curMaster = -1;
                curMasterProc = -1;
                curSlave = -1;
                curSlaveProc = -1;
                curSlaveProcDonorId = -1;

                // Note.  Signalling in global clustering requires
                // me to recognise clustering from separate
                // processors as separate.  In the first phase,
                // this will be used to recognise cluster from
                // each processor as separate and in the second
                // phase it will be used to filter local processor
                // faces from the global patch.  Currently, I am
                // calculating unique cluster index as:
                //
                // id = cluster + procOffset*myProcID
                // HJ, 1/Apr/2009

                // Get live cell restrict
                if (fringeFaceFlips[ffI])
                {
                    // Face pointing into live cell: ie neighbour is live
                    curMaster = localRestrictAddressing[nei[curFace]];
                }
                else
                {
                    // Face pointing out of live cell
                    curMaster = localRestrictAddressing[own[curFace]];
                }

                // Donor cell data is accessible from fringe
                // fringeFaceCells gives acceptor index
                curMasterProc = Pstream::myProcNo();

                curSlave = neighbourRestrictAddressing[fringeFaceCells[ffI]];

                curSlaveProc = acceptorCellDonorProc[fringeFaceCells[ffI]];

                curSlaveProcDonorId =
                    acceptorCellDonorIndex[fringeFaceCells[ffI]];

                // Pout<< "A: " << curMaster << tab << curMasterProc << tab << curSlave << tab << curSlaveProc << tab << curSlaveProcDonorId << endl;
                // Code in current master and slave
                curMaster += procOffset*curMasterProc;
                curSlave += procOffset*curSlaveProc;

                // Look for the master cell.  If it has already got a face,
                // add the coefficient to the face.  If not, create a new
                // face
                if (neighboursTable.found(curMaster))
                {
                    // This master side face already exists

                    // Check all current neighbours to see if the current
                    // slave already exists.  If so, add the coefficient.

                    DynamicList<long, 4>& curNbrs =
                        neighboursTable.find(curMaster)();

                    DynamicList<label, 4>& curNbrsProc =
                        nbrsProcTable.find(curMaster)();

                    DynamicList<label, 4>& curNbrsProcDonorId =
                        nbrsProcDonorIdTable.find(curMaster)();

                    DynamicList<DynamicList<label, 4>, 4>& curFaceFaces =
                        faceFaceTable.find(curMaster)();

                    // Search for coded neighbour
                    bool nbrFound = false;

                    forAll (curNbrs, curNbrI)
                    {
                        // Check neighbour slave
                        if (curNbrs[curNbrI] == curSlave)
                        {
                            nbrFound = true;
                            curFaceFaces[curNbrI].append(ffI);

                            // New agglomeration pair found in already
                            // existing pair
                            nAgglomPairs++;

                            break;
                        }
                    }

                    if (!nbrFound)
                    {
                        curNbrs.append(curSlave);
                        curNbrsProc.append(curSlaveProc);
                        curNbrsProcDonorId.append(curSlaveProcDonorId);

                        DynamicList<label, 4> newFF;
                        newFF.append(ffI);
                        curFaceFaces.append(newFF);

                        // New coarse face created for an existing master
                        nCoarseFaces++;
                        nAgglomPairs++;
                    }
                }
                else
                {
                    // This master has got no neighbours yet.
                    // Add a neighbour, proc and a coefficient as a
                    // new list, thus creating a new face
                    DynamicList<long, 4> newNbrs;
                    newNbrs.append(curSlave);
                    neighboursTable.insert
                    (
                        curMaster,
                        newNbrs
                    );

                    DynamicList<label, 4> newNbrsProc;
                    newNbrsProc.append(curSlaveProc);
                    nbrsProcTable.insert
                    (
                        curMaster,
                        newNbrsProc
                    );

                    DynamicList<label, 4> newNbrsProcDonorId;
                    newNbrsProcDonorId.append(curSlaveProcDonorId);
                    nbrsProcDonorIdTable.insert
                    (
                        curMaster,
                        newNbrsProcDonorId
                    );

                    DynamicList<DynamicList<label, 4>, 4> newFF;
                    newFF.append(DynamicList<label, 4>());
                    newFF[0].append(ffI);
                    faceFaceTable.insert
                    (
                        curMaster,
                        newFF
                    );

                    // New coarse face created for a new master
                    nCoarseFaces++;
                    nAgglomPairs++;
                }
            } // end for all current neighbours
        } // end for all fine faces
    }
    else
    {
        // Coarse level, addressing is stored in faceCells

        // Perform analysis only for local faces
        // HJ, 22/Jun/2016

        label curMasterProc, curSlaveProc, curSlaveProcDonorId;
        long curMaster, curSlave;

        // Size check.  Remove on debug
        if
        (
            localRestrictAddressing.size()
         != fineOversetInterface_.interfaceSize()
        )
        {
            FatalErrorInFunction
                << "Size check failed.  localRestrictAddressing: "
                << localRestrictAddressing.size()
                << " faceCells: "
                << fineOversetInterface_.interfaceSize()
                << abort(FatalError);
        }

        // On the coarse level, acceptor cells are in faceCells_
        // Donor data is served in the same order as acceptors: this is a
        // matched interface
        const labelList& fineAcceptorCells =
            fineOversetInterface_.acceptorCells();

        forAll (fineAcceptorCells, ffI)
        {
            curMaster = -1;
            curMasterProc = -1;
            curSlave = -1;
            curSlaveProc = -1;
            curSlaveProcDonorId = -1;

            // Note.  Signalling in global clustering requires
            // me to recognise clustering from separate
            // processors as separate.  In the first phase,
            // this will be used to recognise cluster from
            // each processor as separate and in the second
            // phase it will be used to filter local processor
            // faces from the global patch.  Currently, I am
            // calculating unique cluster index as:
            //
            // id = cluster + procOffset*myProcID
            // HJ, 1/Apr/2009

            // Not sure if master-slave switching is needed here:
            // each side keeps acceptors as local faceCells and they
            // communicate not in pars but on a more complex pattern

            // Master side
            curMaster = localRestrictAddressing[fineAcceptorCells[ffI]];
            curMasterProc = Pstream::myProcNo();
            curSlave = neighbourRestrictAddressing[ffI];
            curSlaveProc = acceptorCellDonorProc[ffI];
            curSlaveProcDonorId = acceptorCellDonorIndex[ffI];

            // Pout<< "B: " << curMaster << tab << curMasterProc << tab << curSlave << tab << curSlaveProc << tab << curSlaveProcDonorId << endl;
            // Code in current master and slave
            curMaster += procOffset*curMasterProc;
            curSlave += procOffset*curSlaveProc;

            // Look for the master cell.  If it has already got a face,
            // add the coefficient to the face.  If not, create a new face.
            if (neighboursTable.found(curMaster))
            {
                // This master side face already exists

                // Check all current neighbours to see if the current slave
                // already exists and if so, add the fine face
                // to the agglomeration.

                DynamicList<long, 4>& curNbrs =
                    neighboursTable.find(curMaster)();

                DynamicList<label, 4>& curNbrsProc =
                    nbrsProcTable.find(curMaster)();

                DynamicList<label, 4>& curNbrsProcDonorId =
                    nbrsProcDonorIdTable.find(curMaster)();

                DynamicList<DynamicList<label, 4>, 4>& curFaceFaces =
                    faceFaceTable.find(curMaster)();

                // Search for coded neighbour
                bool nbrFound = false;

                forAll (curNbrs, curNbrI)
                {
                    // Check neighbour slave
                    if (curNbrs[curNbrI] == curSlave)
                    {
                        nbrFound = true;
                        curFaceFaces[curNbrI].append(ffI);

                        // New agglomeration pair found in already
                        // existing pair
                        nAgglomPairs++;

                        break;
                    }
                }

                if (!nbrFound)
                {
                    curNbrs.append(curSlave);
                    curNbrsProc.append(curSlaveProc);
                    curNbrsProcDonorId.append(curSlaveProcDonorId);

                    DynamicList<label, 4> newFF;
                    newFF.append(ffI);
                    curFaceFaces.append(newFF);

                    // New coarse face created for an existing master
                    nCoarseFaces++;
                    nAgglomPairs++;
                }
            }
            else
            {
                // This master has got no neighbours yet.  Add a neighbour
                // and a coefficient, thus creating a new face
                DynamicList<long, 4> newNbrs;
                newNbrs.append(curSlave);
                neighboursTable.insert
                (
                    curMaster,
                    newNbrs
                );

                DynamicList<label, 4> newNbrsProc;
                newNbrsProc.append(curSlaveProc);
                nbrsProcTable.insert
                (
                    curMaster,
                    newNbrsProc
                );

                DynamicList<label, 4> newNbrsProcDonorId;
                newNbrsProcDonorId.append(curSlaveProcDonorId);
                nbrsProcDonorIdTable.insert
                (
                    curMaster,
                    newNbrsProcDonorId
                );

                DynamicList<DynamicList<label, 4>, 4> newFF;
                newFF.append(DynamicList<label, 4>());
                newFF[0].append(ffI);
                faceFaceTable.insert
                (
                    curMaster,
                    newFF
                );

                // New coarse face created for a new master
                nCoarseFaces++;
                nAgglomPairs++;
            }
        } // end for all fine faces
    } // end of else in fine level (coarse level)


    // Resize the lists
    faceCells_.setSize(nCoarseFaces);
    fineAddressing_.setSize(nAgglomPairs, -1);
    restrictAddressing_.setSize(nAgglomPairs, -11);

    // All weights are equal to 1: integral matching
    restrictWeights_.setSize(nAgglomPairs, 1.0);

    // Record donor processor for each acceptor (in faceCells_)
    // to establish communication
    donorCellForAcceptor_.setSize(nCoarseFaces);
    donorProcForAcceptor_.setSize(nCoarseFaces);
    donorProcIndexForAcceptor_.setSize(nCoarseFaces);

    List<long> contents = neighboursTable.sortedToc();

    // Note:
    // Each patch has both donors and acceptors.
    // Acceptors are master side: faceCells_ needs to address into the
    // acceptors

    // Donors are on the slave side: this is the data that will need to be
    // sent accross to the donor

    // Reset face counter for re-use
    nCoarseFaces = 0;
    nAgglomPairs = 0;

    // On master side, the owner addressing is stored in table of contents
    forAll (contents, masterI)
    {
        const label curMaster = contents[masterI];

        const DynamicList<long, 4>& curNbrs =
            neighboursTable.find(curMaster)();

        const DynamicList<label, 4>& curNbrsProc =
            nbrsProcTable.find(curMaster)();

        const DynamicList<label, 4>& curNbrsProcDonorId =
            nbrsProcDonorIdTable.find(curMaster)();

        const DynamicList<DynamicList<label, 4>, 4>& curFaceFaces =
            faceFaceTable.find(curMaster)();

        forAll (curNbrs, curNbrI)
        {
            // Record new coarse acceptor
            faceCells_[nCoarseFaces] =
                curMaster - procOffset*Pstream::myProcNo();

            // Record donor for coarse acceptor (can be on different processor)
            donorCellForAcceptor_[nCoarseFaces] = curNbrs[curNbrI];

            // Record donor processor
            donorProcForAcceptor_[nCoarseFaces] = curNbrsProc[curNbrI];

            // Record donor processir index for acceptor
            donorProcIndexForAcceptor_[nCoarseFaces] =
                curNbrsProcDonorId[curNbrI];

            // Get faces and weights
            const DynamicList<label, 4>& facesIter = curFaceFaces[curNbrI];

            forAll (facesIter, facesIterI)
            {
                fineAddressing_[nAgglomPairs] = facesIter[facesIterI];
                restrictAddressing_[nAgglomPairs] = nCoarseFaces;

                nAgglomPairs++;
            }

            nCoarseFaces++;
        }
    }

    // Info<< "Finished level.  Sizes " << nl
    //     << "faceCells_: " << faceCells_
    //     << " fineAddressing_: " << fineAddressing_.size()
    //     << " restrictAddressing_: " << restrictAddressing_.size()
    //     << " donorCells_: " << donorCells_.size()
    //     << " donorCellForAcceptor_: " << donorCellForAcceptor_.size()
    //     << " donorProcIndexForAcceptor_: " << donorProcIndexForAcceptor_.size()
    //     << endl;

    if (Pstream::parRun())
    {
        initMap();
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::oversetAMGInterface::~oversetAMGInterface()
{
    deleteDemandDrivenData(mapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::oversetAMGInterface::interfaceInternalField
(
    const unallocLabelList& iF
) const
{
    // Return complete field: needed for both donors and acceptors
    // HJ, 16/Sep/2019
    return tmp<labelField>(new labelField(iF));
}


Foam::tmp<Foam::scalarField> Foam::oversetAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    // Escape agglomerating internal coefficients: zero size
    // HJ, 16/Sep/2019
    if (fineCoeffs.empty())
    {
        return tmp<scalarField>(new scalarField());
    }

    tmp<scalarField> tcoarseCoeffs(new scalarField(size(), 0.0));
    scalarField& coarseCoeffs = tcoarseCoeffs();

    // Added weights to account for non-integral matching
    forAll (restrictAddressing_, ffi)
    {
        coarseCoeffs[restrictAddressing_[ffi]] +=
            restrictWeights_[ffi]*fineCoeffs[fineAddressing_[ffi]];
    }

    return tcoarseCoeffs;
}


bool Foam::oversetAMGInterface::master() const
{
    return fineOversetInterface_.master();
}


bool Foam::oversetAMGInterface::fineLevel() const
{
    return false;
}


Foam::label Foam::oversetAMGInterface::interfaceSize() const
{
    return interfaceSize_;
}


const Foam::lduMesh& Foam::oversetAMGInterface::ldu() const
{
    return AMGInterface::ldu();
}


const Foam::oversetMesh& Foam::oversetAMGInterface::overset() const
{
    // Overset should not be accessed from coarse levels
    FatalErrorInFunction
        << "Requested fine addressing at coarse level"
            << abort(FatalError);

    // Dummy return
    return fineOversetInterface_.overset();
}


const Foam::mapDistribute& Foam::oversetAMGInterface::map() const
{
    if (!mapPtr_)
    {
        initMap();
    }

    return *mapPtr_;
}


void Foam::oversetAMGInterface::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{}


Foam::tmp<Foam::labelField> Foam::oversetAMGInterface::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    return labelField::null();
}


void Foam::oversetAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    // Repackage donor data to acceptors
    // Note: result is copied as internal field, re-scaled and passed across
    if (Pstream::parRun())
    {
        labelTransferBuffer_ = iF;

        map().distribute(labelTransferBuffer_);
    }
    else
    {
        labelTransferBuffer_ = labelField(iF, donorCellForAcceptor_);
    }
}


Foam::tmp<Foam::labelField> Foam::oversetAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return tmp<labelField>(new labelField(labelTransferBuffer_));
}


void Foam::oversetAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const scalarField& iF
) const
{
    // Repackage donor data to acceptors
    // Note: result is copied as internal field, re-scaled and passed across

    // Repackage donor data to acceptors
    // Note: result is copied as internal field, re-scaled and passed across
    if (Pstream::parRun())
    {
        fieldTransferBuffer_ = iF;

        map().distribute(fieldTransferBuffer_);
    }
    else
    {
        fieldTransferBuffer_ = scalarField(iF, donorCellForAcceptor_);
    }
}


Foam::tmp<Foam::scalarField> Foam::oversetAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const scalarField&
) const
{
    return tmp<scalarField>(new scalarField(fieldTransferBuffer_));
}


// ************************************************************************* //
