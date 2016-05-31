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
    if (!Pstream::parRun())
    {
        FatalErrorIn("void ggiAMGInterface::initFastReduce() const")
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // Init should be executed only once
    initReduce_ = true;

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

    reduce(zoneProcID, maxOp<labelField>());

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

    labelListList shadowToReceiveAddr(Pstream::nProcs());

    // Get the list of what my shadow needs to receive from my zone
    // on all other processors
    shadowToReceiveAddr[Pstream::myProcNo()] = shadowZa;
    Pstream::gatherList(shadowToReceiveAddr);
    Pstream::scatterList(shadowToReceiveAddr);

    // Now local zone indices contain the index of a local face that will
    // provide the data.  For faces that are not local, the index will be -1
 
    // Find out where my zone data is going to

    // Make a sending sub-map
    // It tells me which data I will send to which processor
    labelListList sendMap(Pstream::nProcs());

    // Collect local labels to be sent to each processor
    forAll (shadowToReceiveAddr, procI)
    {
        const labelList& curProcSend = shadowToReceiveAddr[procI];

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

    // Map will return the object of the size of remote zone
    // HJ, 9/May/2016
    mapPtr_ = new mapDistribute(zoneSize(), sendMap, constructMap);
}


const Foam::mapDistribute& Foam::ggiAMGInterface::map() const
{
    if (!mapPtr_)
    {
        initFastReduce();
    }

    return *mapPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiAMGInterface::ggiAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    AMGInterface(lduMesh),
    fineGgiInterface_(refCast<const ggiLduInterface>(fineInterface)),
    zoneSize_(0),
    zoneAddressing_(),
    initReduce_(false),
    mapPtr_(NULL)
{
    // Note.
    // All processors will do the same coarsening and then filter
    // the addressing to the local processor
    // HJ, 1/Apr/2009

    // Note.
    // Signalling in global clustering requires me to recognise clustering
    // from separate processors as separate.  In the first phase, this will be
    // used to recognise cluster from each processor as separate and in the
    // second phase it will be used to filter local processor faces from
    // the global patch.
    // Currently, I am calculating unique cluster index as:
    //
    // id = cluster + procOffset*myProcID
    // With procOffset = 1 million, this should be sufficient for 2000 CPUs
    // with 2 million coarse cells each.  For larger numbers, I need a
    // larger max int, which can be changed on request
    // HJ, 1/Apr/2009

    // Expand the local and neighbour addressing to full zone size
    labelField localExpandAddressing(fineGgiInterface_.zoneSize(), 0);

    // Memory management, local
    {
        const labelList& addr = fineGgiInterface_.zoneAddressing();

        forAll (addr, i)
        {
            localExpandAddressing[addr[i]] =
                localRestrictAddressing[i] + procOffset*Pstream::myProcNo();
        }

        if (!localParallel())
        {
            reduce(localExpandAddressing, sumOp<labelField>());
        }
    }

    labelField neighbourExpandAddressing
    (
        fineGgiInterface_.shadowInterface().zoneSize(),
        0
    );

    // Memory management, neighbour
    {
        const labelList& addr =
            fineGgiInterface_.shadowInterface().zoneAddressing();

        forAll (addr, i)
        {
            neighbourExpandAddressing[addr[i]] =
                neighbourRestrictAddressing[i]
              + procOffset*Pstream::myProcNo();
        }

        if (!localParallel())
        {
            reduce(neighbourExpandAddressing, sumOp<labelField>());
        }
    }

    // Make a lookup table of entries for owner/neighbour.
    // All sizes are guessed at the size of fine interface
    // HJ, 19/Feb/2009

    HashTable<SLList<label>, label, Hash<label> > neighboursTable
    (
        localExpandAddressing.size()
    );

    // Table of face-sets to be agglomerated
    HashTable<SLList<SLList<label> >, label, Hash<label> > faceFaceTable
    (
        localExpandAddressing.size()
    );

    // Table of face-sets weights to be agglomerated
    HashTable<SLList<SLList<scalar> >, label, Hash<label> >
        faceFaceWeightsTable
        (
            localExpandAddressing.size()
        );

    // Count the number of coarse faces
    label nCoarseFaces = 0;

    // Count the number of agglomeration pairs
    label nAgglomPairs = 0;

    // On the fine level, addressing is made in a labelListList
    if (fineGgiInterface_.fineLevel())
    {
        const labelListList& fineAddr = fineGgiInterface_.addressing();
        const scalarListList& fineWeights = fineGgiInterface_.weights();

        forAll (fineAddr, ffI)
        {
            const labelList& curFineNbrs = fineAddr[ffI];
            const scalarList& curFineWeigts = fineWeights[ffI];

            forAll (curFineNbrs, nbrI)
            {
                label curMaster = -1;
                label curSlave = -1;

                // My label = ffI
                // Nbr label = nnI
                const label nnI = curFineNbrs[nbrI];
                const scalar curNW = curFineWeigts[nbrI];

                if (fineGgiInterface_.master())
                {
                    // Master side
                    curMaster = localExpandAddressing[ffI];
                    curSlave = neighbourExpandAddressing[nnI];
                }
                else
                {
                    // Slave side
                    curMaster = neighbourExpandAddressing[nnI];
                    curSlave = localExpandAddressing[ffI];
                }

                // Look for the master cell.  If it has already got a face,
                // add the coefficient to the face.  If not, create a new
                // face.
                if (neighboursTable.found(curMaster))
                {
                    // Check all current neighbours to see if the current
                    // slave already exists.  If so, add the coefficient.

                    SLList<label>& curNbrs = neighboursTable.find(curMaster)();

                    SLList<SLList<label> >& curFaceFaces =
                        faceFaceTable.find(curMaster)();

                    SLList<SLList<scalar> >& curFaceWeights =
                        faceFaceWeightsTable.find(curMaster)();

                    bool nbrFound = false;

                    SLList<label>::iterator nbrsIter = curNbrs.begin();

                    SLList<SLList<label> >::iterator faceFacesIter =
                        curFaceFaces.begin();

                    SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                        curFaceWeights.begin();

                    for
                    (
                        ;
                        nbrsIter != curNbrs.end()
                     && faceFacesIter != curFaceFaces.end()
                     && faceFaceWeightsIter != curFaceWeights.end();
                        ++nbrsIter, ++faceFacesIter, ++faceFaceWeightsIter
                    )
                    {
                        // Check neighbour slave
                        if (nbrsIter() == curSlave)
                        {
                            nbrFound = true;
                            faceFacesIter().append(ffI);
                            faceFaceWeightsIter().append(curNW);
                            nAgglomPairs++;

                            break;
                        }
                    }

                    if (!nbrFound)
                    {
                        curNbrs.append(curSlave);
                        curFaceFaces.append(SLList<label>(ffI));
                        curFaceWeights.append(SLList<scalar>(curNW));

                        // New coarse face created
                        nCoarseFaces++;
                        nAgglomPairs++;
                    }
                }
                else
                {
                    // This master has got no neighbours yet.  Add a neighbour
                    // and a coefficient, thus creating a new face
                    neighboursTable.insert(curMaster, SLList<label>(curSlave));

                    faceFaceTable.insert
                    (
                        curMaster,
                        SLList<SLList<label> >(SLList<label>(ffI))
                    );

                    faceFaceWeightsTable.insert
                    (
                        curMaster,
                        SLList<SLList<scalar> >(SLList<scalar>(curNW))
                    );

                    // New coarse face created
                    nCoarseFaces++;
                    nAgglomPairs++;
                }
            } // end for all current neighbours
        } // end for all fine faces
    }
    else
    {
        // Coarse level, addressing is stored in faceCells
        forAll (localExpandAddressing, ffi)
        {
            label curMaster = -1;
            label curSlave = -1;

            // Do switching on master/slave indices based on the
            // owner/neighbour of the processor index such that
            // both sides get the same answer.
            if (master())
            {
                // Master side
                curMaster = localExpandAddressing[ffi];
                curSlave = neighbourExpandAddressing[ffi];
            }
            else
            {
                // Slave side
                curMaster = neighbourExpandAddressing[ffi];
                curSlave = localExpandAddressing[ffi];
            }

            // Look for the master cell.  If it has already got a face,
            // add the coefficient to the face.  If not, create a new face.
            if (neighboursTable.found(curMaster))
            {
                // Check all current neighbours to see if the current slave
                // already exists and if so, add the fine face
                // to the agglomeration.

                SLList<label>& curNbrs = neighboursTable.find(curMaster)();

                SLList<SLList<label> >& curFaceFaces =
                    faceFaceTable.find(curMaster)();

                    SLList<SLList<scalar> >& curFaceWeights =
                        faceFaceWeightsTable.find(curMaster)();

                bool nbrFound = false;

                SLList<label>::iterator nbrsIter = curNbrs.begin();

                SLList<SLList<label> >::iterator faceFacesIter =
                    curFaceFaces.begin();

                SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                    curFaceWeights.begin();

                for
                (
                    ;
                    nbrsIter != curNbrs.end()
                 && faceFacesIter != curFaceFaces.end()
                 && faceFaceWeightsIter != curFaceWeights.end();
                    ++nbrsIter, ++faceFacesIter, ++faceFaceWeightsIter
                )
                {
                    // Check neighbour slave
                    if (nbrsIter() == curSlave)
                    {
                        nbrFound = true;
                        faceFacesIter().append(ffi);
                        // Add dummy weight
                        faceFaceWeightsIter().append(1.0);
                        nAgglomPairs++;
                        break;
                    }
                }

                if (!nbrFound)
                {
                    curNbrs.append(curSlave);
                    curFaceFaces.append(SLList<label>(ffi));
                    // Add dummy weight
                    curFaceWeights.append(SLList<scalar>(1.0));

                    // New coarse face created
                    nCoarseFaces++;
                    nAgglomPairs++;
                }
            }
            else
            {
                // This master has got no neighbours yet.  Add a neighbour
                // and a coefficient, thus creating a new face
                neighboursTable.insert(curMaster, SLList<label>(curSlave));

                faceFaceTable.insert
                (
                    curMaster,
                    SLList<SLList<label> >(SLList<label>(ffi))
                );

                // Add dummy weight
                faceFaceWeightsTable.insert
                (
                    curMaster,
                    SLList<SLList<scalar> >(SLList<scalar>(1.0))
                );

                // New coarse face created
                nCoarseFaces++;
                nAgglomPairs++;
            }
        } // end for all fine faces
    }

    faceCells_.setSize(nCoarseFaces, -1);
    fineAddressing_.setSize(nAgglomPairs, -1);
    restrictAddressing_.setSize(nAgglomPairs, -1);
    restrictWeights_.setSize(nAgglomPairs);

    labelList contents = neighboursTable.toc();

    // Sort makes sure the order is identical on both sides.
    // HJ, 20/Feb.2009
    sort(contents);

    // Grab zone size and create zone addressing
    zoneSize_ = nCoarseFaces;

    zoneAddressing_.setSize(nCoarseFaces);
    label nProcFaces = 0;

    // Reset face counter for re-use
    nCoarseFaces = 0;
    nAgglomPairs = 0;

    if (master())
    {
        // On master side, the owner addressing is stored in table of contents
        forAll (contents, masterI)
        {
            SLList<label>& curNbrs =
                neighboursTable.find(contents[masterI])();

            // Note: neighbour processor index is irrelevant.  HJ, 1/Apr/2009

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            SLList<SLList<scalar> >& curFaceWeights =
                faceFaceWeightsTable.find(contents[masterI])();

            SLList<label>::iterator nbrsIter = curNbrs.begin();
            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                curFaceWeights.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end()
             && faceFacesIter != curFaceFaces.end()
             && faceFaceWeightsIter != curFaceWeights.end();
                ++nbrsIter, ++faceFacesIter, ++faceFaceWeightsIter
            )
            {
                // Check if master is on local processor
                if
                (
                    contents[masterI] >= procOffset*Pstream::myProcNo()
                 && contents[masterI] < procOffset*(Pstream::myProcNo() + 1)
                )
                {
                    // Record that this face belongs locally
                    zoneAddressing_[nProcFaces] = nCoarseFaces;
                    faceCells_[nProcFaces] =
                        contents[masterI] - procOffset*Pstream::myProcNo();
                    nProcFaces++;

                    SLList<label>::iterator facesIter =
                        faceFacesIter().begin();
                    SLList<scalar>::iterator weightsIter =
                        faceFaceWeightsIter().begin();

                    for
                    (
                        ;
                        facesIter != faceFacesIter().end()
                     && weightsIter != faceFaceWeightsIter().end();
                        ++facesIter, ++weightsIter
                    )
                    {
                        fineAddressing_[nAgglomPairs] = facesIter();
                        restrictAddressing_[nAgglomPairs] = nCoarseFaces;
                        restrictWeights_[nAgglomPairs] = weightsIter();
                        nAgglomPairs++;
                    }
                }

                // Not a local face, but still created in global zone
                nCoarseFaces++;
            }
        }

        // Resize arrays: not all of ggi is used locally
        faceCells_.setSize(nProcFaces);
        zoneAddressing_.setSize(nProcFaces);

        fineAddressing_.setSize(nAgglomPairs);
        restrictAddressing_.setSize(nAgglomPairs);
        restrictWeights_.setSize(nAgglomPairs);
    }
    else
    {
        // On slave side, the owner addressing is stored in linked lists
        forAll (contents, masterI)
        {
            // Note: master processor index is irrelevant.  HJ, 1/Apr/2009

            SLList<label>& curNbrs = neighboursTable.find(contents[masterI])();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            SLList<SLList<scalar> >& curFaceWeights =
                faceFaceWeightsTable.find(contents[masterI])();

            SLList<label>::iterator nbrsIter = curNbrs.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            SLList<SLList<scalar> >::iterator faceFaceWeightsIter =
                curFaceWeights.begin();
            for
            (
                ;
                nbrsIter != curNbrs.end()
             && faceFacesIter != curFaceFaces.end()
             && faceFaceWeightsIter != curFaceWeights.end();
                ++nbrsIter, ++faceFacesIter, ++faceFaceWeightsIter
            )
            {
                // Check if the face is on local processor
                if
                (
                    nbrsIter() >= procOffset*Pstream::myProcNo()
                 && nbrsIter() < procOffset*(Pstream::myProcNo() + 1)
                )
                {
                    // Record that this face belongs locally.
                    zoneAddressing_[nProcFaces] = nCoarseFaces;
                    faceCells_[nProcFaces] =
                        nbrsIter() - procOffset*Pstream::myProcNo();
                    nProcFaces++;

                    SLList<label>::iterator facesIter =
                        faceFacesIter().begin();

                    SLList<scalar>::iterator weightsIter =
                        faceFaceWeightsIter().begin();

                    for
                    (
                        ;
                        facesIter != faceFacesIter().end()
                     && weightsIter != faceFaceWeightsIter().end();
                        ++facesIter, ++weightsIter
                    )
                    {
                        fineAddressing_[nAgglomPairs] = facesIter();
                        restrictAddressing_[nAgglomPairs] = nCoarseFaces;
                        restrictWeights_[nAgglomPairs] = weightsIter();
                        nAgglomPairs++;
                    }
                }

                nCoarseFaces++;
            }
        }

        // Resize arrays: not all of ggi is used locally
        faceCells_.setSize(nProcFaces);
        zoneAddressing_.setSize(nProcFaces);
        fineAddressing_.setSize(nAgglomPairs);
        restrictAddressing_.setSize(nAgglomPairs);
        restrictWeights_.setSize(nAgglomPairs);
    }
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
    // Note: reconsider better parallel communication here.
    // Currently expanding to full zone size
    // HJ, 16/Mar/2016

    // Reassemble fine coefficients to full fine zone size
    scalarField zoneFineCoeffs(fineGgiInterface_.zoneSize(), 0);

    const labelList& fineZa = fineGgiInterface_.zoneAddressing();

    forAll (fineZa, i)
    {
        zoneFineCoeffs[fineZa[i]] = fineCoeffs[i];
    }

    // Reduce zone data
    if (!localParallel())
    {
        reduce(zoneFineCoeffs, sumOp<scalarField>());
    }

    scalarField zoneCoarseCoeffs(zoneSize(), 0);

    // Restrict coefficient
    forAll(restrictAddressing_, ffi)
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


// ************************************************************************* //
