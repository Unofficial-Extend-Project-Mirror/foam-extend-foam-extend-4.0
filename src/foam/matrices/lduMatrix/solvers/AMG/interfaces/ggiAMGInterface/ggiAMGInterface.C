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

\*---------------------------------------------------------------------------*/

#include "ggiAMGInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        AMGInterface,
        ggiAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ggiAMGInterface::initFastReduce() const
{
    if (mapPtr_)
    {
        FatalErrorIn("void ggiAMGInterface::initFastReduce() const")
            << "map already calculated"
            << abort(FatalError);
    }

    if (!Pstream::parRun())
    {
        FatalErrorIn("void ggiAMGInterface::initFastReduce() const")
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // If the processor is not in the GGI comm, create a dummy map
    if (Pstream::myProcNo(comm()) == -1)
    {
        mapPtr_ = new mapDistribute
        (
            zoneSize(),  // This is zero if there are no local GGI faces
            labelListList(Pstream::nProcs()),
            labelListList(Pstream::nProcs())
        );

        return;
    }
    Info<< "Start initFastReduce " << lTime_.elapsedCpuTime() << endl;
    // From here on, work on processors within the communicator
    // HJ, 20/Sep/2016

    // Note: this is different from ggiPolyPatch comms because zone
    // is the same on master the slave side.
    // HJ, 31/May/2016

    // Establish parallel comms pattern

    // Get addressing
    const labelList& za = zoneAddressing();
    const labelList& shadowZa = shadowInterface().zoneAddressing();

    // Make a zone-sized field and fill it in with proc markings for processor
    // that holds and requires the data
    labelField zoneProcID(zoneSize(), -1);

    forAll (za, zaI)
    {
        zoneProcID[za[zaI]] = Pstream::myProcNo();
    }

    // Note: reduce with a comm will only be present on processors containing
    // master or slave faces.  Other processors created a dummy map above
    // HJ, 20/Sep/2016
    reduce(zoneProcID, maxOp<labelField>(), tag(), comm());

    // Find out where my zone data is coming from
    labelList nRecv(Pstream::nProcs(), 0);

    forAll (shadowZa, shadowZaI)
    {
        nRecv[zoneProcID[shadowZa[shadowZaI]]]++;
    }

    // Make a receiving sub-map
    // It tells me which data I will receive from which processor and
    // where I need to put it into the remoteZone data before the mapping
    labelListList constructMap(Pstream::nProcs());

    // Size the receiving list
    forAll (nRecv, procI)
    {
        constructMap[procI].setSize(nRecv[procI]);
    }

    // Reset counters for processors
    nRecv = 0;

    forAll (shadowZa, shadowZaI)
    {
        label recvProc = zoneProcID[shadowZa[shadowZaI]];

        constructMap[recvProc][nRecv[recvProc]] = shadowZa[shadowZaI];

        nRecv[recvProc]++;
    }

    // Make the sending sub-map
    // It tells me which data is required from me to be sent to which
    // processor

    // Algorithm
    // - expand the local zone faces with indices into a size of local zone
    // - go through remote zone addressing on all processors
    // - find out who hits my faces
    labelList localZoneIndices(zoneSize(), -1);

    forAll (za, zaI)
    {
        localZoneIndices[za[zaI]] = zaI;
    }

    // Note: rewrite this to use comm().  However, all proc indicators will be
    // done with local comms addressing instead of the global processor number
    // HJ, 20/Sep/2016
    labelListList shadowToReceiveAddr(Pstream::nProcs(comm()));

    // Get the list of what my shadow needs to receive from my zone
    // on all other processors
    shadowToReceiveAddr[Pstream::myProcNo(comm())] = shadowZa;

    // Note gather-scatter with a comm.  Processors without GGI faces
    // have been eliminated beforehand
    Pstream::gatherList(shadowToReceiveAddr, tag(), comm());
    Pstream::scatterList(shadowToReceiveAddr, tag(), comm());

    // Now local zone indices contain the index of a local face that will
    // provide the data.  For faces that are not local, the index will be -1

    // Find out where my zone data is going to

    // Make a sending sub-map
    // It tells me which data I will send to which processor
    labelListList sendMap(Pstream::nProcs());

    // Collect local labels to be sent to each processor
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        // Look for processors that are in the current comm, ie. containing
        // faces from the ggi patch
        if (Pstream::procNo(comm(), procI) != -1)
        {
            const labelList& curProcSend =
                shadowToReceiveAddr[Pstream::procNo(comm(), procI)];

            // Find out how much of my data is going to this processor
            label nProcSend = 0;

            forAll (curProcSend, sendI)
            {
                if (localZoneIndices[curProcSend[sendI]] > -1)
                {
                    nProcSend++;
                }
            }

            if (nProcSend > 0)
            {
                // Collect the indices
                labelList& curSendMap = sendMap[procI];

                curSendMap.setSize(nProcSend);

                // Reset counter
                nProcSend = 0;

                forAll (curProcSend, sendI)
                {
                    if (localZoneIndices[curProcSend[sendI]] > -1)
                    {
                        curSendMap[nProcSend] =
                            localZoneIndices[curProcSend[sendI]];
                        nProcSend++;
                    }
                }
            }
        }
    }

    // Map will return the object of the size of remote zone
    // HJ, 9/May/2016
    mapPtr_ = new mapDistribute(zoneSize(), sendMap, constructMap);
    Info<< "End initFastReduce " << lTime_.elapsedCpuTime() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiAMGInterface::ggiAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    AMGInterface(lduMesh),
    fineGgiInterface_(refCast<const ggiLduInterface>(fineInterface)),
    zoneSize_(0),
    zoneAddressing_(),
    procMasterFaces_(),
    comm_(fineGgiInterface_.comm()),
    tag_(fineGgiInterface_.tag()),
    mapPtr_(NULL),
    lTime_()
{
    // New algorithm will assemble local clusters on the master side and
    // create zone ordering by collecting all faces (coarse pairs) from proc0,
    // followed by proc 1 etc.  This avoids global communication and allows
    // each processor only to perform the analysis on locally created coarse
    // faces
    // HJ, 13/Jun/2016
    Info<< "Start ggiAMGInterface constructor for size "
        << fineGgiInterface_.interfaceSize() << ": " 
        << lTime_.elapsedCpuTime() << endl;

    // Note: local addressing contains only local faces
    const labelList& fineZa =  fineGgiInterface_.zoneAddressing();

    // Create addressing for neighbour faces.  Note: expandAddrToZone will
    // expand the addressing to zone size, including communications.
    // Faces which are not used locally will be marked by NaNs
    // HJ, 13/Jun/2016
    labelField neighbourExpandAddressing = neighbourRestrictAddressing;

    // Expand neighbour side to get all the data required from other
    // processors.
    // Note: neigbour is now the size of remote zone
    if (!fineGgiInterface_.shadowInterface().localParallel())
    {
        fineGgiInterface_.shadowInterface().expandAddrToZone
        (
            neighbourExpandAddressing
        );
    }

    // Create addressing for neighbour processors.  Note: expandAddrToZone will
    // expand the addressing to zone size.  HJ, 13/Jun/2016
    labelField neighbourExpandProc
    (
        fineGgiInterface_.shadowInterface().interfaceSize(),
        Pstream::myProcNo()
    );

    // Expand neighbour side to get all the data required from other
    // processors.
    // Note: neigbour is now the size of remote zone
    if (!fineGgiInterface_.shadowInterface().localParallel())
    {
        fineGgiInterface_.shadowInterface().expandAddrToZone
        (
            neighbourExpandProc
        );
    }

    // Note: neighbourExpandAddressing and neighbourExpandProc
    // will be filled with NaNs for faces which are not local

    // Make a lookup table of entries for owner/neighbour.
    // All sizes are guessed at the size of fine interface
    // HJ, 19/Feb/2009

    // Note: Guessing size of HashTable to fine interface size

    // Coded neighbour index. Note: using long int to simplify encoding
    // HJ, 1/Aug/2016
    HashTable<SLList<long>, long, Hash<long> > neighboursTable
    (
        Foam::max(128, fineGgiInterface_.interfaceSize()/4)
    );

    // Neignbour processor index
    HashTable<SLList<label>, label, Hash<label> > nbrsProcTable
    (
        Foam::max(128, fineGgiInterface_.interfaceSize()/4)
    );

    // Neighbour face-faces addressing for a face with split neighbours
    HashTable<SLList<SLList<label> >, label, Hash<label> > faceFaceTable
    (
        Foam::max(128, fineGgiInterface_.interfaceSize()/4)
    );

    HashTable<SLList<SLList<label> >, label, Hash<label> > faceFaceNbrTable
    (
        Foam::max(128, fineGgiInterface_.interfaceSize()/4)
    );

    // Neighbour face-faces weights for a face with split neighbours
    HashTable<SLList<SLList<scalar> >, label, Hash<label> >
        faceFaceWeightsTable
        (
            Foam::max(128, fineGgiInterface_.interfaceSize()/4)
        );

    // Count the number of coarse faces
    label nCoarseFaces = 0;

    // Count the number of agglomeration pairs
    label nAgglomPairs = 0;

    // On the fine level, addressing is made in a labelListList
    if (fineGgiInterface_.fineLevel())
    {
        // This addressing defines how to interpolate for all zone faces
        // across the interface
        const labelListList& fineAddr = fineGgiInterface_.addressing();
        const scalarListList& fineWeights = fineGgiInterface_.weights();

        // Perform analysis only for local faces
        // HJ, 22/Jun/2016
        forAll (fineZa, fineZaI)
        {
            // Get the local face (from zone) to analyse
            const label ffI = fineZa[fineZaI];

            const labelList& curFineNbrs = fineAddr[ffI];
            const scalarList& curFineWeigts = fineWeights[ffI];

            forAll (curFineNbrs, nbrI)
            {
                long curMaster = -1;
                label curMasterProc = -1;
                long curSlave = -1;
                label curSlaveProc = -1;

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
                //
                // With procOffset = 1 million, this should be
                // sufficient for 2000 CPUs with 2 million coarse
                // cells each.  For larger numbers, I need a
                // larger max int, which can be changed on request
                // HJ, 1/Apr/2009

                // My label = ffI
                // Nbr label = nnI
                const label nnI = curFineNbrs[nbrI];
                const scalar curNW = curFineWeigts[nbrI];

                if (fineGgiInterface_.master())
                {
                    // Master side
                    curMaster = localRestrictAddressing[fineZaI];
                    curMasterProc = Pstream::myProcNo();
                    curSlave = neighbourExpandAddressing[nnI];
                    curSlaveProc = neighbourExpandProc[nnI];
                }
                else
                {
                    curMaster = neighbourExpandAddressing[nnI];
                    curMasterProc = neighbourExpandProc[nnI];
                    curSlave = localRestrictAddressing[fineZaI];
                    curSlaveProc = Pstream::myProcNo();
                }

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

                    SLList<long>& curNbrs =
                        neighboursTable.find(curMaster)();

                    SLList<label>& curNbrsProc =
                        nbrsProcTable.find(curMaster)();

                    SLList<SLList<label> >& curFaceFaces =
                        faceFaceTable.find(curMaster)();

                    SLList<SLList<label> >& curFaceFaceNbrs =
                        faceFaceNbrTable.find(curMaster)();

                    SLList<SLList<scalar> >& curFaceWeights =
                        faceFaceWeightsTable.find(curMaster)();

                    // Search for coded neighbour
                    bool nbrFound = false;

                    SLList<long>::iterator nbrsIter = curNbrs.begin();

                    SLList<label>::iterator nbrsProcIter =
                        curNbrsProc.begin();

                    SLList<SLList<label> >::iterator faceFacesIter =
                        curFaceFaces.begin();

                    SLList<SLList<label> >::iterator faceFaceNbrsIter =
                        curFaceFaceNbrs.begin();

                    SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                    curFaceWeights.begin();

                    for
                    (
                        ;
                        nbrsIter != curNbrs.end()
                     && nbrsProcIter != curNbrsProc.end()
                     && faceFacesIter != curFaceFaces.end()
                     && faceFaceNbrsIter != curFaceFaceNbrs.end()
                     && faceFaceWeightsIter != curFaceWeights.end();
                        ++nbrsIter, ++nbrsProcIter,
                        ++faceFacesIter, ++faceFaceNbrsIter,
                        ++faceFaceWeightsIter
                    )
                    {
                        // Check neighbour slave
                        if (nbrsIter() == curSlave)
                        {
                            nbrFound = true;
                            faceFacesIter().append(ffI);
                            faceFaceNbrsIter().append(nbrI);
                            faceFaceWeightsIter().append(curNW);

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
                        curFaceFaces.append(SLList<label>(ffI));
                        curFaceFaceNbrs.append(SLList<label>(nbrI));
                        curFaceWeights.append(SLList<scalar>(curNW));

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
                    neighboursTable.insert
                    (
                        curMaster,
                        SLList<long>(curSlave)
                    );

                    nbrsProcTable.insert
                    (
                        curMaster,
                        SLList<label>(curSlaveProc)
                    );

                    faceFaceTable.insert
                    (
                        curMaster,
                        SLList<SLList<label> >(SLList<label>(ffI))
                    );

                    faceFaceNbrTable.insert
                    (
                        curMaster,
                        SLList<SLList<label> >(SLList<label>(nbrI))
                    );

                    faceFaceWeightsTable.insert
                    (
                        curMaster,
                        SLList<SLList<scalar> >(SLList<scalar>(curNW))
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
        // This addressing defines which faces from zone are local

        // Perform analysis only for local faces
        // HJ, 22/Jun/2016
        forAll (fineZa, fineZaI)
        {
            // Get the local face (from zone) to analyse
            const label ffI = fineZa[fineZaI];

            long curMaster = -1;
            label curMasterProc = -1;
            long curSlave = -1;
            label curSlaveProc = -1;

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
            //
            // With procOffset = 1 million, this should be
            // sufficient for 2000 CPUs with 2 million coarse
            // cells each.  For larger numbers, I need a
            // larger max int, which can be changed on request
            // HJ, 1/Apr/2009

            if (fineGgiInterface_.master())
            {
                // Master side
                curMaster = localRestrictAddressing[fineZaI];
                curMasterProc = Pstream::myProcNo();
                curSlave = neighbourExpandAddressing[ffI];
                curSlaveProc = neighbourExpandProc[ffI];
            }
            else
            {
                curMaster = neighbourExpandAddressing[ffI];
                curMasterProc = neighbourExpandProc[ffI];
                curSlave = localRestrictAddressing[fineZaI];
                curSlaveProc = Pstream::myProcNo();
            }

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

                SLList<long>& curNbrs = neighboursTable.find(curMaster)();

                SLList<label>& curNbrsProc =
                    nbrsProcTable.find(curMaster)();

                SLList<SLList<label> >& curFaceFaces =
                    faceFaceTable.find(curMaster)();

                SLList<SLList<label> >& curFaceFaceNbrs =
                    faceFaceNbrTable.find(curMaster)();

                SLList<SLList<scalar> >& curFaceWeights =
                    faceFaceWeightsTable.find(curMaster)();

                bool nbrFound = false;

                SLList<long>::iterator nbrsIter = curNbrs.begin();

                SLList<label>::iterator nbrsProcIter =
                    curNbrsProc.begin();

                SLList<SLList<label> >::iterator faceFacesIter =
                    curFaceFaces.begin();

                SLList<SLList<label> >::iterator faceFaceNbrsIter =
                    curFaceFaceNbrs.begin();

                SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                    curFaceWeights.begin();

                for
                (
                    ;
                    nbrsIter != curNbrs.end()
                 && faceFacesIter != curFaceFaces.end()
                 && faceFaceNbrsIter != curFaceFaceNbrs.end()
                 && faceFaceWeightsIter != curFaceWeights.end();
                    ++nbrsIter, ++nbrsProcIter,
                    ++faceFacesIter, ++faceFaceNbrsIter, ++faceFaceWeightsIter
                )
                {
                    // Check neighbour slave
                    if (nbrsIter() == curSlave)
                    {
                        nbrFound = true;
                        faceFacesIter().append(ffI);
                        // Add dummy nbr
                        faceFaceNbrsIter().append(0);
                        // Add dummy weight
                        faceFaceWeightsIter().append(1.0);

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
                    curFaceFaces.append(SLList<label>(ffI));
                    // Add dummy nbr
                    curFaceFaceNbrs.append(SLList<label>(0));
                    // Add dummy weight
                    curFaceWeights.append(SLList<scalar>(1.0));

                    // New coarse face created for an existing master
                    nCoarseFaces++;
                    nAgglomPairs++;
                }
            }
            else
            {
                // This master has got no neighbours yet.  Add a neighbour
                // and a coefficient, thus creating a new face
                neighboursTable.insert
                (
                    curMaster,
                    SLList<long>(curSlave)
                );

                nbrsProcTable.insert
                (
                    curMaster,
                    SLList<label>(curSlaveProc)
                );

                faceFaceTable.insert
                (
                    curMaster,
                    SLList<SLList<label> >(SLList<label>(ffI))
                );

                // Add dummy nbr
                faceFaceNbrTable.insert
                (
                    curMaster,
                    SLList<SLList<label> >(SLList<label>(0))
                );

                // Add dummy weight
                faceFaceWeightsTable.insert
                (
                    curMaster,
                    SLList<SLList<scalar> >(SLList<scalar>(1.0))
                );

                // New coarse face created for a new master
                nCoarseFaces++;
                nAgglomPairs++;
            }
        } // end for all fine faces
    } // end of else in fine level (coarse level)

    // Since only local faces are analysed, lists can now be resized
    faceCells_.setSize(nCoarseFaces, -1);
    fineAddressing_.setSize(nAgglomPairs, -1);
    restrictAddressing_.setSize(nAgglomPairs, -1);
    restrictWeights_.setSize(nAgglomPairs);

    // In order to assemble the coarse global face zone, find out
    // how many faces have been created on each processor.
    // Note that masters and slaves both count faces so we will
    //  only ask master sizes to count
    labelList nCoarseFacesPerProc(Pstream::nProcs(), 0);

    nCoarseFacesPerProc[Pstream::myProcNo()] = nCoarseFaces;

    // Reduce number of coarse faces per proc
    // Note: reducing with comms means that the processors with no
    // contact with the GGI interface will have zero zone size.
    // This needs to be handled separately in the initFastReduce
    // HJ, 20/Sep/2016
    reduce(nCoarseFacesPerProc, sumOp<labelList>(), tag(), comm());

    // Coarse global face zone is assembled by adding all faces from proc0,
    // followed by all faces from proc1 etc.
    // Therefore, on procN, my master offset
    // will be equal to the sum of numbers of coarse faces on all
    // processors before mine
    // HJ, 13/Jun/2016

    label coarseGlobalFaceOffset = 0;

    for (label i = 0; i < Pstream::myProcNo(); i++)
    {
        coarseGlobalFaceOffset += nCoarseFacesPerProc[i];
    }

    // Grab zone size and create zone addressing
    zoneSize_ = sum(nCoarseFacesPerProc);

    zoneAddressing_.setSize(nCoarseFaces);

    // Both master and slave have done agglomeration, but only master
    // will construct the global faces index.
    // To avoid searching, master will prepare a list of global face
    // indices that appear on each processor in order and communicate them
    // to the slave
    // The slave will then know which is the next global face from which
    // processor and pick them out in the same order

    // Global faces shall be assembled by the increasing label of master
    // cluster ID.
    List<long> contents = neighboursTable.toc();

    // Sort makes sure the order is identical on both sides.
    // HJ, 20/Feb/2009 and 6/Jun/2016
    sort(contents);

    // Note: Restriction is done on master side only because this is where
    // the local zone is created.  HJ, 1/Aug/2016
    if (master())
    {
        // Note:
        // When I am agglomerating the master, faces are stacked up in order
        // but on the slave side, all I know is the master cluster index and
        // not a master coarse face index.  Therefore:
        // - master needs to be agglomerated first
        // - once master is agglomerated, I need to signal to the slave side
        //   the global coarse face zone index

        // For each new global face created on master proc,
        // record its index under the slave proc array
        List<SLList<label> > procMasterFacesLL(Pstream::nProcs());

        // Note: zone addressing will be assembled only for local clusters
        // using the coarseGlobalFaceOffset
        // HJ, 13/Jun/2016
        label nProcFaces = 0;

        // Reset face counter for re-use
        nCoarseFaces = 0;
        nAgglomPairs = 0;

        // Note:
        // Since clustering has now happened only on local faces,
        // addressing and all other array work on local indices and
        // not on the coarse global zone
        // HJ, 13/Jun/2016

        // On master side, the owner addressing is stored in table of contents
        forAll (contents, masterI)
        {
            SLList<long>& curNbrs =
                neighboursTable.find(contents[masterI])();

            SLList<label>& curNbrsProc =
                nbrsProcTable.find(contents[masterI])();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            SLList<SLList<scalar> >& curFaceWeights =
                faceFaceWeightsTable.find(contents[masterI])();

            SLList<long>::iterator nbrsIter = curNbrs.begin();

            SLList<label>::iterator nbrsProcIter =
                curNbrsProc.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                curFaceWeights.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end()
             && nbrsProcIter != curNbrsProc.end()
             && faceFacesIter != curFaceFaces.end()
             && faceFaceWeightsIter != curFaceWeights.end();
                ++nbrsIter, ++nbrsProcIter,
                ++faceFacesIter, ++faceFaceWeightsIter
            )
            {
                // Check if master is on local processor: no longer needed,
                // as only local processor is being searched.  HJ, 13/Jun/2016

                // Get faces and weights
                SLList<label>::iterator facesIter =
                    faceFacesIter().begin();

                SLList<scalar>::iterator weightsIter =
                    faceFaceWeightsIter().begin();

                // Record that this face belongs locally
                // Use offset to indicate its position in the list
                zoneAddressing_[nProcFaces] =
                    nProcFaces + coarseGlobalFaceOffset;

                // Record master cluster index
                faceCells_[nProcFaces] =
                    contents[masterI] - procOffset*Pstream::myProcNo();

                // Record global processor face
                procMasterFacesLL[nbrsProcIter()].append
                (
                    nProcFaces + coarseGlobalFaceOffset
                );

                // Collect agglomeration data
                for
                (
                    ;
                    facesIter != faceFacesIter().end()
                 && weightsIter != faceFaceWeightsIter().end();
                    ++facesIter, ++weightsIter
                )
                {
                    fineAddressing_[nAgglomPairs] = facesIter();

                    // Master processor zone face is calculated from
                    // global offset
                    restrictAddressing_[nAgglomPairs] =
                        nProcFaces + coarseGlobalFaceOffset;

                    restrictWeights_[nAgglomPairs] = weightsIter();
                    nAgglomPairs++;
                }

                nProcFaces++;
            }
        }

        // No need to resize arrays only local faces are used
        // HJ, 1/Aug/2016

        // Re-pack singly linked list of processor master faces
        // and pass to other processors

        procMasterFaces_.setSize(Pstream::nProcs());

        // Copy self
        procMasterFaces_[Pstream::myProcNo()] =
            procMasterFacesLL[Pstream::myProcNo()];

        if (Pstream::parRun())
        {
            const List<labelPair> schedule =
                fineGgiInterface_.map().schedule();

            // Do the comms
            forAll (schedule, i)
            {
                const label sendProc = schedule[i].first();
                const label recvProc = schedule[i].second();

                if (Pstream::myProcNo() == sendProc)
                {
                    OPstream toNbr
                    (
                        Pstream::scheduled,
                        Pstream::procNo(comm(), recvProc),
                        0,
                        tag(),
                        comm()
                    );
                    toNbr << labelList(procMasterFacesLL[recvProc]);
                }
                else if (Pstream::myProcNo() == recvProc)
                {
                    IPstream fromNbr
                    (
                        Pstream::scheduled,
                        Pstream::procNo(comm(), sendProc),
                        0,
                        tag(),
                        comm()
                    );

                    procMasterFaces_[sendProc] = labelList(fromNbr);
                }
                else
                {
                    FatalErrorIn("...")
                        << "My proc number " << Pstream::myProcNo()
                            << " is neither a sender nor a receiver: "
                            << schedule[i]
                            << abort(FatalError);
                }
            }
        }
        Info<< "ggiAMGInterface end agglom master " << lTime_.elapsedCpuTime() << endl;
    }
    // Agglomerate slave
    else
    {
        // Note: zone addressing will be assembled only for local clusters
        // using the coarseGlobalFaceOffset
        // HJ, 13/Jun/2016
        label nProcFaces = 0;

        // Get master side procMasterFaces
        // Note: this needs to be picked up from coarse interfaces rather than
        // the matrix, as the coarse matrix assembly is not complete yet.
        // Further, master has got a lower index in the list, meaning that
        // it has already completed the assembly
        const ggiAMGInterface& shadowGGI =
            refCast<const ggiAMGInterface>(coarseInterfaces[shadowIndex()]);

        const labelListList& masterProcMasterFaces =
            shadowGGI.procMasterFaces();

        // Reset face counter for re-use
        nCoarseFaces = 0;
        nAgglomPairs = 0;

        // Count how many global faces are used for each proc on the other side
        labelList npmf(Pstream::nProcs(), 0);

        // On slave side, the owner addressing is stored in linked lists
        forAll (contents, masterI)
        {
            SLList<long>& curNbrs = neighboursTable.find(contents[masterI])();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            SLList<SLList<label> >& curFaceFaceNbrs =
                faceFaceNbrTable.find(contents[masterI])();

            SLList<SLList<scalar> >& curFaceWeights =
                faceFaceWeightsTable.find(contents[masterI])();

            SLList<long>::iterator nbrsIter = curNbrs.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            SLList<SLList<label> >::iterator faceFaceNbrsIter =
                curFaceFaceNbrs.begin();

            SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                curFaceWeights.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end()
             && faceFacesIter != curFaceFaces.end()
             && faceFaceNbrsIter != curFaceFaceNbrs.end()
             && faceFaceWeightsIter != curFaceWeights.end();
                ++nbrsIter, ++faceFacesIter, ++faceFaceNbrsIter,
                ++faceFaceWeightsIter
            )
            {
                // Check if the face is on local processor: no longer needed,
                // as only local processor is being searched.  HJ, 13/Jun/2016

                SLList<label>::iterator facesIter =
                    faceFacesIter().begin();

                SLList<label>::iterator faceNbrsIter =
                    faceFaceNbrsIter().begin();

                SLList<scalar>::iterator weightsIter =
                    faceFaceWeightsIter().begin();

                // Find neighbour proc index from the first face
                // on the other side
                label nbrProc;

                if (fineGgiInterface_.fineLevel())
                {
                    const labelListList& fineAddr =
                        fineGgiInterface_.addressing();

                    nbrProc = neighbourExpandProc
                        [fineAddr[facesIter()][faceNbrsIter()]];
                }
                else
                {
                    nbrProc = neighbourExpandProc[facesIter()];
                }

                // Read coarse face index
                const label coarseFace =
                    masterProcMasterFaces[nbrProc][npmf[nbrProc]];

                // and mark it as used
                npmf[nbrProc]++;

                zoneAddressing_[nProcFaces] = coarseFace;
                faceCells_[nProcFaces] =
                    nbrsIter() - procOffset*Pstream::myProcNo();

                for
                (
                    ;
                    facesIter != faceFacesIter().end()
                 && weightsIter != faceFaceWeightsIter().end();
                    ++facesIter, ++weightsIter
                )
                {
                    fineAddressing_[nAgglomPairs] = facesIter();
                    restrictAddressing_[nAgglomPairs] = coarseFace;
                    restrictWeights_[nAgglomPairs] = weightsIter();

                    nAgglomPairs++;
                }

                nProcFaces++;
            }
        }
        Info<< "ggiAMGInterface end agglom slave " << lTime_.elapsedCpuTime() << endl;
    }
    Info<< "End ggiAMGInterface constructor " << lTime_.elapsedCpuTime() << endl;
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::ggiAMGInterface::~ggiAMGInterface()
{
    deleteDemandDrivenData(mapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::ggiAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    // HJ, HERE THIS SHOULD BE REMOVED: NO LONGER NEEDED BECAUSE ALL ADDRESSING IS LOCAL

    // Note: reconsider better parallel communication here.
    // Currently expanding to full zone size
    // HJ, 16/Mar/2016

    // Reassemble fine coefficients to full fine zone size
    // No need to initialise to zero, as only local coefficients
    // are used.  HJ, 9/Jun/2016
    scalarField zoneFineCoeffs(fineGgiInterface_.zoneSize());

    const labelList& fineZa = fineGgiInterface_.zoneAddressing();

    forAll (fineZa, i)
    {
        zoneFineCoeffs[fineZa[i]] = fineCoeffs[i];
    }

    // Reduce zone data is not required: all coefficients are local
    // HJ, 9/Jun/2016

    scalarField zoneCoarseCoeffs(zoneSize(), 0);

    // Restrict coefficient
    forAll (restrictAddressing_, ffi)
    {
        zoneCoarseCoeffs[restrictAddressing_[ffi]] +=
            restrictWeights_[ffi]*zoneFineCoeffs[fineAddressing_[ffi]];
    }

    tmp<scalarField> tcoarseCoeffs(new scalarField(size(), 0.0));
    scalarField& coarseCoeffs = tcoarseCoeffs();

    // Filter zone coefficients to local field
    const labelList& za = zoneAddressing();

    forAll (za, i)
    {
        coarseCoeffs[i] = zoneCoarseCoeffs[za[i]];
    }

    return tcoarseCoeffs;
}


bool Foam::ggiAMGInterface::master() const
{
    return fineGgiInterface_.master();
}


bool Foam::ggiAMGInterface::fineLevel() const
{
    return false;
}


Foam::label Foam::ggiAMGInterface::shadowIndex() const
{
    return fineGgiInterface_.shadowIndex();
}


const Foam::ggiLduInterface& Foam::ggiAMGInterface::shadowInterface() const
{
    return refCast<const ggiLduInterface>(ldu().interfaces()[shadowIndex()]);
}


Foam::label Foam::ggiAMGInterface::interfaceSize() const
{
    return faceCells_.size();
}

Foam::label Foam::ggiAMGInterface::zoneSize() const
{
    return zoneSize_;
}


const Foam::labelList& Foam::ggiAMGInterface::zoneAddressing() const
{
    return zoneAddressing_;
}


const Foam::labelListList& Foam::ggiAMGInterface::addressing() const
{
    FatalErrorIn("const labelListList& ggiAMGInterface::addressing() const")
        << "Requested fine addressing at coarse level"
        << abort(FatalError);

    return labelListList::null();
}


bool Foam::ggiAMGInterface::localParallel() const
{
    return fineGgiInterface_.localParallel();
}


const Foam::labelListList& Foam::ggiAMGInterface::procMasterFaces() const
{
    if (!master())
    {
        FatalErrorIn
        (
            "const labelListList& ggiGAMGInterface::procMasterFaces() const"
        )   << "Requester procMasterFaces from a slave.  This is not allowed"
            << abort(FatalError);
    }

    return procMasterFaces_;
}


const Foam::mapDistribute& Foam::ggiAMGInterface::map() const
{
    if (!mapPtr_)
    {
        initFastReduce();
    }

    return *mapPtr_;
}


const Foam::scalarListList& Foam::ggiAMGInterface::weights() const
{
    FatalErrorIn("const labelListList& ggiAMGInterface::weights() const")
        << "Requested fine addressing at coarse level"
        << abort(FatalError);

    return scalarListList::null();
}


const Foam::tensorField& Foam::ggiAMGInterface::forwardT() const
{
    return fineGgiInterface_.forwardT();
}


const Foam::tensorField& Foam::ggiAMGInterface::reverseT() const
{
    return fineGgiInterface_.reverseT();
}


void Foam::ggiAMGInterface::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    // Label transfer is local
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::ggiAMGInterface::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    // Label transfer is local without global reduction
    return this->shadowInterface().labelTransferBuffer();
}


void Foam::ggiAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    // Label transfer is local without global reduction
    labelTransferBuffer_ = interfaceInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::ggiAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList&
) const
{
    return shadowInterface().labelTransferBuffer();
}


void Foam::ggiAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const scalarField& iF
) const
{
    scalarField pif = interfaceInternalField(iF);

    // New treatment.  HJ, 26/Jun/2011
    fieldTransferBuffer_ = fastReduce(pif);
}


Foam::tmp<Foam::scalarField> Foam::ggiAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const scalarField&
) const
{
    return shadowInterface().fieldTransferBuffer();
}


void Foam::ggiAMGInterface::expandAddrToZone(labelField& lf) const
{
    if (!localParallel())
    {
        lf = fastExpand(lf);
    }
}


// ************************************************************************* //
