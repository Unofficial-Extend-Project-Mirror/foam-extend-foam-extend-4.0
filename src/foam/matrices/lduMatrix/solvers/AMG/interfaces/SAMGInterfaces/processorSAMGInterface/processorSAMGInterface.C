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

#include "processorSAMGInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorSAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        SAMGInterface,
        processorSAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorSAMGInterface::processorSAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const crMatrix& prolongation,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const crMatrix& nbrInterfaceProlongation
)
:
    SAMGInterface(lduMesh, prolongation, nbrInterfaceProlongation),
    fineProcInterface_(refCast<const processorLduInterface>(fineInterface)),
    comm_(fineProcInterface_.comm()),
    tag_(fineProcInterface_.tag())
{
    /*
    // On master:
    //  - receive slave's filtered prolongation
    //  - create a transpose of your local prolongation
    //  - do the triple product: store the coefficients and the addresses
    //    row = row of restriction, column = column of prolongation
    //  - the coefficients should be stored row by row, but with randomized
    //    columns
    //  - send the addressing to slave to sort its coefficients to match
    //    the ordering on master (is it possible to sort at the same time?)

    // Needed:
    //  - labelList masterOwner,
    //  - labelList masterNeighbour,
    //  - label nCoarseCoeffs


    // Question: Should we put the masterOwner and masterNeighbour into a linked
    // list which can be searched quickly? In this case we can create the
    // addressing on master, send and at the same time, calculate the coarse
    // coefficients on both sides and store into corresponding order

    // Question: How can I save the addressing of the boundary coeffs for
    // the next level?

    // Question: With such boundary coefficients, my vector-matrix product is a
    // matrix-matrix product - change multiplication rule in SAMGInterfaceFields
    // or account for it with addressing?
//------------------------------------------------------------------------------

    // Count coarse coefficients and create addressing on master
    if (myProcNo() < neighbProcNo())
    {
//------------------------------------------------------------------------------
//                              COUNT COARSE COEFFS
//------------------------------------------------------------------------------
        // Filtered prolongation matrix from my side
        const labelList& localOwner = prolongation.crAddr().column();

        // Filtered prolongation matrix from my side
        tmp<crMatrix> tProlongationT(prolongation.T());
        const labelList& rRowStart = tProlongationT->crAddr().rowStart();
        const labelList& rColumn = tProlongationT->crAddr().column();
        const label rNRows = tProlongationT->crAddr().nRows();

        // Filtered prolongation matrix from neighbour - to obtain restriction,
        // make a transpose
        const labelList& pRowStart =
            nbrInterfaceProlongation.crAddr().rowStart();
        const labelList& pColumn = nbrInterfaceProlongation.crAddr().column();
        const label pNCols = nbrInterfaceProlongation.crAddr().nCols();

        labelList coeffMark(pNCols, -1);
        label nCoarseCoeffs = 0;

        // Row of R
        for (label ir = 0; ir < rNRows; ir++)
        {
            //Restart coeffMark for each restriction row
            coeffMark = -1;

            // Row start addressing R
            for
            (
                label indexR = rRowStart[ir];
                indexR < rRowStart[ir + 1];
                indexR++
            )
            {
                // Column of coefficient in R
                const label jr = rColumn[indexR];

                for
                (
                    label indexP = pRowStart[jr];
                    indexP < pRowStart[jr + 1];
                    indexP++
                )
                {
                    label jp = pColumn[indexP];

                    if (coeffMark[jp] == -1)
                    {
                        // Found a new coarse coeff!
                        coeffMark[jp] = nCoarseCoeffs;
                        nCoarseCoeffs++;
                    }
                }
            }
        }

//------------------------------------------------------------------------------
//                      CREATE HASH TABLE FOR ADDRESSING
//------------------------------------------------------------------------------

        HashTable<label, labelPair, Hash<labelPair> > coarseLevel(nCoarseCoeffs);
        faceCells_.setSize(nCoarseCoeffs);

        // Reset for creating addressing
        nCoarseCoeffs = 0;
        coeffMark = -1;

        // Row of R
        for (label ir = 0; ir < rNRows; ir++)
        {
            // Restart coeffMark for each restriction row
            coeffMark = -1;

            // Row start addressing R
            for
            (
                label indexR = rRowStart[ir];
                indexR < rRowStart[ir + 1];
                indexR++
            )
            {
                // Column of coefficient in R
                const label jr = rColumn[indexR];

                // Grab column of the original prolongation, this is the
                // faceCell_!
                const label ja = localOwner[ir];

                for
                (
                    label indexP = pRowStart[jr];
                    indexP < pRowStart[jr + 1];
                    indexP++
                )
                {
                    label jp = pColumn[indexP];

                    label identify = coeffMark[jp];

                    if (identify == -1)
                    {
                        // Found a new coarse coeff!
                        identify = nCoarseCoeffs;
                        coeffMark[jp] = nCoarseCoeffs;

                        faceCells_[identify] = ja;

                        // Save address and coeff into HashTable
                        coarseLevel.insert(labelPair(ir, jp), nCoarseCoeffs);
                        nCoarseCoeffs++;
                    }
                }
            }
        }

        // Send to slave
        OPstream toNbr(Pstream::blocking, neighbProcNo());

        toNbr<< coarseLevel;

    }
    // Slave
    else
    {
        // Slave must receive the hash table sent by the master to know how to
        // sort its boundary coefficients

        IPstream fromNbr(Pstream::blocking, neighbProcNo());

        masterCoarseLevel_ =
            HashTable<label, labelPair, Hash<labelPair> >(fromNbr);

        // Filtered prolongation matrix from my side
        const labelList& localOwner = prolongation.crAddr().column();

        // Filtered prolongation matrix from my side
        tmp<crMatrix> tProlongationT(prolongation.T());
        const labelList& rRowStart = tProlongationT->crAddr().rowStart();
        const labelList& rColumn = tProlongationT->crAddr().column();
        const label rNRows = tProlongationT->crAddr().nRows();

        // Filtered prolongation matrix from neighbour - to obtain restriction,
        // make
        // a transpose
        const labelList& pRowStart = nbrInterfaceProlongation.crAddr().rowStart();
        const labelList& pColumn = nbrInterfaceProlongation.crAddr().column();
        const label pNCols = nbrInterfaceProlongation.crAddr().nCols();

        faceCells_.setSize(masterCoarseLevel_.size());

        labelList coeffMark(pNCols, -1);

        // Row of R
        for (label ir = 0; ir < rNRows; ir++)
        {
            // Restart coeffMark for each restriction row
            coeffMark = -1;

            // Row start addressing R
            for
            (
                label indexR = rRowStart[ir];
                indexR < rRowStart[ir + 1];
                indexR++
            )
            {
                // Column of coefficient in R
                const label jr = rColumn[indexR];

                // Grab column of the original prolongation -  this is the
                // faceCell_!!!
                label ja = localOwner[ir];

                for
                (
                    label indexP = pRowStart[jr];
                    indexP < pRowStart[jr + 1];
                    indexP++
                )
                {
                    label jp = pColumn[indexP];

                    label identify = coeffMark[jp];

                    if (identify == -1)
                    {
                        label address = masterCoarseLevel_[labelPair(jp,ja)];

                        // Found a new coarse coeff!
                        identify = address;

                        coeffMark[jp] = address;
                        faceCells_[identify] = ja;
                    }
                }
            }
        }
    }


    // Analyse the local and neighbour row label:
    //  local coarse, remote coarse = regular coarse face
    //  local coarse, remote fine = local expanded face: receive prolonged data
    //  local fine, remote coarse = neighbour expanded face: send prolonged data
    //  local fine, remote fine = local and neighbour expanded face: send and
    //                            receive

    // Algorithm
    // Go through local row labels and examine cases
    // 1) coarse local and coarse remote:
    //    create coarse processor face and set faceCells.  Set weight to 1
    // 2) coarse local and fine remote:
    //    will receive coarse neighbours and weights from opposite side
    // 3) fine local and coarse remote:
    //    assemble local coarse neighbours and weights from the prolongation
    // 4) fine local and fine remote - ignore
    //
    // On completion of selection:
    //    send and receive local coarse neighbours for fine local/coarse remote
    //    sort coarse equations by increasing coarse index from master side
    //    create faceCells, fineAddressing and fineWeights

    // First collect and communicate internal neighbours from the local fine
    // side (ie neighbour processor cell is coarse)
    HashTable<labelList, label, Hash<label> > neighboursFromLocalFine
    (
        Foam::max(128, fineProcInterface_.interfaceSize()/4)
    );

    HashTable<scalarField, label, Hash<label> > weightsFromLocalFine
    (
        Foam::max(128, fineProcInterface_.interfaceSize()/4)
    );

    // For FINE to FINE communication - local boundary coefficient multiplyed by
    // local weight - send to the other side for triple product
    HashTable<labelList, label, Hash<label> > neighboursFromLocalFineToFine
    (
        Foam::max(128, fineProcInterface_.interfaceSize()/4)
    );

    HashTable<scalarField, label, Hash<label> > weightsFromLocalFineToFine
    (
        Foam::max(128, fineProcInterface_.interfaceSize()/4)
    );


    // Get access to the prolongation addressing and coefficients
    const labelList& pRowStart = prolongation.crAddr().rowStart();
    const labelList& pColumn = prolongation.crAddr().column();
    const scalarField& pCoeffs = prolongation.coeffs();

    // Get fine faceCells
    const labelList& fineFaceCells = fineInterface.faceCells();

    // Collect local fine to neighbour coarse connections for communication
    // TU, 19 May 2017: and fine to fine (communication through triple product)
    forAll (localRowLabel, faceI)
    {
        if
        (
            localRowLabel[faceI] < 0
         && neighbourRowLabel[faceI] >= 0
        )
        {
            // Found local fine to neighbour coarse interface

            // Collect local coarse neighbours and weights from the
            // prolongation matrix to send to other processor
            const label curStart = pRowStart[fineFaceCells[faceI]];
            const label curEnd = pRowStart[fineFaceCells[faceI] + 1];
            const label nCoarse = curEnd - curStart;

            labelList nbrs(nCoarse);
            scalarField weights(nCoarse);

            forAll (nbrs, i)
            {
                nbrs[i] = pColumn[curStart + i];
                weights[i] = pCoeffs[curStart + i];
            }

            // Insert neighbours under remote coarse index
            neighboursFromLocalFine.insert(neighbourRowLabel[faceI], nbrs);
            weightsFromLocalFine.insert(neighbourRowLabel[faceI], weights);
        }
        // FINE to FINE communication - TU, May 2017
        if
        (
            localRowLabel[faceI] < 0
         && neighbourRowLabel[faceI] < 0
        )
        {
            // Found local FINE to neighbour FINE interface

            // First, multiply all local interface boundary coefficients with
            // local prolongation weights - you will send this to the other side
            // in localBoundaryProlongation

            // Collect local coarse neighbours and weights from the
            // prolongation matrix to send to other processor
            const label curStart = pRowStart[fineFaceCells[faceI]];
            const label curEnd = pRowStart[fineFaceCells[faceI] + 1];
            const label nCoeffs = curEnd - curStart;

            labelList nbrs(nCoeffs);
            scalarField weights(nCoeffs);

            forAll (nbrs, i)
            {
                nbrs[i] = pColumn[curStart + i];
                weights[i] = pCoeffs[curStart + i];
            }

            // Insert neighbours under remote coarse index
            neighboursFromLocalFineToFine.insert(neighbourRowLabel[faceI], nbrs);
            weightsFromLocalFineToFine.insert(neighbourRowLabel[faceI], weights);
        }
    }

    // Receive remote prolongation data
    HashTable<labelList, label, Hash<label> > neighboursFromRemoteFine;
    HashTable<scalarField, label, Hash<label> > weightsFromRemoteFine;
    HashTable<labelList, label, Hash<label> > neighboursFromRemoteFineToFine;
    HashTable<scalarField, label, Hash<label> > weightsFromRemoteFineToFine;

    // Send and receive the addressing from the other side
    {
        OPstream toNbr(Pstream::blocking, neighbProcNo());
        toNbr<< neighboursFromLocalFine << weightsFromLocalFine;
    }

    {
        IPstream fromNbr(Pstream::blocking, neighbProcNo());

        neighboursFromRemoteFine =
            HashTable<labelList, label, Hash<label> >(fromNbr);

        weightsFromRemoteFine =
            HashTable<scalarField, label, Hash<label> >(fromNbr);

        neighboursFromRemoteFineToFine =
            HashTable<labelList, label, Hash<label> >(fromNbr);

        weightsFromRemoteFineToFine =
            HashTable<scalarField, label, Hash<label> >(fromNbr);

    }

    // Assemble connectivity

    // Resize arrays to size of fine interface
    // Note: it is technically possible to have MORE coarse processor faces
    // so a check will be made in the end
    faceCells_.setSize(5*fineProcInterface_.interfaceSize());
    fineAddressing_.setSize(5*fineProcInterface_.interfaceSize());
    fineWeights_.setSize(5*fineProcInterface_.interfaceSize());

    // Count coarse faces
    label nCoarseFaces = 0;

    // Collect coarse-to-fine connections
    forAll (localRowLabel, faceI)
    {
        if
        (
            localRowLabel[faceI] >= 0
         && neighbourRowLabel[faceI] >= 0
        )
        {
            // Found local coarse to neighbour coarse face

            // Create new coarse face
            faceCells_[nCoarseFaces] = localRowLabel[faceI];
            fineAddressing_[nCoarseFaces] = faceI;
            fineWeights_[nCoarseFaces] = 1;
            nCoarseFaces++;
        }
        else if
        (
            localRowLabel[faceI] < 0
         && neighbourRowLabel[faceI] >= 0
        )
        {
            // Found local fine to neighbour coarse face

            // Pick up local prolongation coarse entries and for all
            // add a new face with appropriate weight
            // Note: faceCells changes due to (internal) local coarse cells
            const labelList& curLocalCoarseNbrs =
                neighboursFromLocalFine[neighbourRowLabel[faceI]];

            const scalarField& curLocalCoarseWeights =
                weightsFromLocalFine[neighbourRowLabel[faceI]];

            forAll (curLocalCoarseNbrs, curNbrI)
            {
                // Create new coarse face
                faceCells_[nCoarseFaces] = curLocalCoarseNbrs[curNbrI];
                fineAddressing_[nCoarseFaces] = faceI;
                fineWeights_[nCoarseFaces] = curLocalCoarseWeights[curNbrI];
                nCoarseFaces++;
            }
        }
        else if
        (
            localRowLabel[faceI] >= 0
         && neighbourRowLabel[faceI] < 0
        )
        {
            // Found local coarse to neighbour fine face

            // Pick up neighbour prolongation coarse entries and for all
            // add a new face with appropriate weight
            // Note: faceCells remains the same: single local coarse cell
            const labelList& curNbrCoarseNbrs =
                neighboursFromRemoteFine[localRowLabel[faceI]];

            const scalarField& curNbrCoarseWeights =
                weightsFromRemoteFine[localRowLabel[faceI]];

            forAll (curNbrCoarseNbrs, curNbrI)
            {
                // Create new coarse face
                faceCells_[nCoarseFaces] = localRowLabel[faceI];
                fineAddressing_[nCoarseFaces] = faceI;
                fineWeights_[nCoarseFaces] = curNbrCoarseWeights[curNbrI];
                nCoarseFaces++;
            }
        }
        else
        {
            // Fine to fine.
            // Establish connection as in triple product: multiply restriction
            // matrix on my side, send to the other side and multiply with
            // prolongation there

            // Get fine faceCells
            const scalarField& fineLocalCoeffs =
            // Pick up local prolongation coarse entries and for all
            // add a new face with appropriate weight
            // Note: faceCells changes, but how???

            const labelList& curLocalCoarseNbrs =
                neighboursFromLocalFineToFine[neighbourRowLabel[faceI]];

            const scalarField& curLocalCoarseWeights =
                weightsFromLocalFineToFine[neighbourRowLabel[faceI]];

            forAll (curLocalCoarseNbrs, curNbrI)
            {
                // Create new coarse face
                faceCells_[nCoarseFaces] = curLocalCoarseNbrs[curNbrI];
                fineAddressing_[nCoarseFaces] = faceI;
                fineWeights_[nCoarseFaces] =
                    curLocalCoarseWeights[curNbrI]*fineLocalCoeffs[faceI];
                nCoarseFaces++;

                // Store fineFaceCells, recording the fine index of prolonged
                // cell locally and all coarse neighbours on the other side
                // HJ, HERE!!!
            }
        }
    }

    // Check for fixed lists
    if (nCoarseFaces > 5*fineProcInterface_.interfaceSize())
    {
        FatalErrorIn("processorSAMGInterface::processorSAMGInterface(...)")
            << "Coarse SAMG processor siginificantly bigger than fine: "
            << "nCoarseFaces = " << nCoarseFaces
            << " nFineFaces = " << fineProcInterface_.interfaceSize()
            << abort(FatalError);
    }

    // Resize arrays to final size
    faceCells_.setSize(nCoarseFaces);
    fineAddressing_.setSize(nCoarseFaces);
    fineWeights_.setSize(nCoarseFaces);
*/
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::processorSAMGInterface::~processorSAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorSAMGInterface::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    send(commsType, interfaceData);
}


Foam::tmp<Foam::labelField> Foam::processorSAMGInterface::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorSAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    send(commsType, interfaceInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorSAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorSAMGInterface::initProlongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    // Send prolongation matrix, using IOstream operators
    OPstream toNbr(Pstream::blocking, neighbProcNo());
    toNbr<< filteredP;
}


Foam::autoPtr<Foam::crMatrix>
Foam::processorSAMGInterface::prolongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    IPstream fromNbr(Pstream::blocking, neighbProcNo());

    autoPtr<crMatrix> tnbrP(new crMatrix(fromNbr));

    return tnbrP;
}


// ************************************************************************* //
