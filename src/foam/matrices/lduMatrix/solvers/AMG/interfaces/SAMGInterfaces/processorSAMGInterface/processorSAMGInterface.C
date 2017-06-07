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
    const crMatrix& interfaceProlongation,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const crMatrix& nbrInterfaceProlongation
)
:
    SAMGInterface(lduMesh, interfaceProlongation, nbrInterfaceProlongation),
    fineProcInterface_(refCast<const processorLduInterface>(fineInterface)),
    comm_(fineProcInterface_.comm()),
    tag_(fineProcInterface_.tag())
{
    // MASTER processor
    if (myProcNo() < neighbProcNo())
    {
        // READ FIRST
//------------------------------------------------------------------------------
        // This code is written for the FILTERED local matrix. Why?
        // Because the addressing is very natural. If this is not ok, filter the
        // big matrix here using fineInterface.faceCells(), transpose, and use
        // the same names for crMatrix arrays (then, there is no need to change
        // the triple product code).
//------------------------------------------------------------------------------

        // Algorithm details:
        // Go through the triple product similar to SAMG Policy - collect
        // coarse addressing on master - send it to slave to sort its addressing
        // accordingly - HashTable containing labelPair(owner, neighbour) as key
        // and coarse face index as entry

        // The resulting contributions from triple product will be saved using
        // the following arrays:
        // - faceCells_: for coarse matrix, it tells us the owner of the
        //               boundary coefficient (on local side)
        // - fineAddressing_: saves the index of the fine boundary coefficient
        //                    in the natural order of appearance on master in
        //                    triple product
        // - restrictAddressing_: saves the coarse face for which the
        //                        contribution (coefficient) needs to sum-up
        // - restrictWeights_: weights for each fine boundary coefficient in
        //                     natural order of appearance on master

        // *Note: (A*B)^T = B^T*A^T, which is exactly what I have. Use the fact
        // on the slave processor!!!

        // Addressing of coarse faces - for sorting faceCells_ on slave
        HashTable<label, labelPair, Hash<labelPair> > coarseLevelTable;

        // Lists for saving addressing and weights - give size, but check at the
        // end!
        DynamicList<label> dynFaceCells(5*fineProcInterface_.interfaceSize());
        DynamicList<label> dynFineAddr(5*fineProcInterface_.interfaceSize());
        DynamicList<scalar> dynRestrictW(5*fineProcInterface_.interfaceSize());
        DynamicList<label> dynRestrictAddr(5*fineProcInterface_.interfaceSize());

        // Filtered prolongation matrix from my side
        crMatrix prolongationT = interfaceProlongation.T();
        const labelList& rRowStart = prolongationT.crAddr().rowStart();
        const labelList& rColumn = prolongationT.crAddr().column();
        const scalarField& rCoeffs = prolongationT.coeffs();
        const label rNRows = prolongationT.crAddr().nRows();

        // Filtered prolongation matrix from neighbour - to obtain restriction,
        // make a transpose
        const labelList& pRowStart =
            nbrInterfaceProlongation.crAddr().rowStart();
        const labelList& pColumn = nbrInterfaceProlongation.crAddr().column();
        const scalarField& pCoeffs = nbrInterfaceProlongation.coeffs();
        const label pNCols = nbrInterfaceProlongation.crAddr().nCols();

        labelList coeffMark(pNCols, -1);
        label nCoarseCoeffs = 0;
        label nCoarseContribs = 0;

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
                // This is an important information: it tells us which
                // boundary coefficient this restriction coeff multiplies! It is
                // the index of the boundary coeff in the boundary array.
                const label jr = rColumn[indexR];

                for
                (
                    label indexP = pRowStart[jr];
                    indexP < pRowStart[jr + 1];
                    indexP++
                )
                {
                    // Column of P, goes into coarse address
                    label jp = pColumn[indexP];

                    // To which coeff do I add myself to? Does the address
                    // in this row already exist?
                    label identify = coeffMark[jp];

                    // If the coeff at this address doesn't exist
                    if (identify == -1)
                    {
                        // Found a new coarse coeff!
                        identify = nCoarseCoeffs;
                        coeffMark[jp] = nCoarseCoeffs;

                        dynFaceCells.append(ir);

                        // Save address and coeff into HashTable for sorting
                        // faceCells on slave
                        coarseLevelTable.insert
                        (
                            labelPair(ir, jp),
                            nCoarseCoeffs
                        );
                        nCoarseCoeffs++;
                    }

                    dynFineAddr.append(jr);
                    dynRestrictW.append(rCoeffs[indexR]*pCoeffs[indexP]);
                    dynRestrictAddr.append(identify);

                    nCoarseContribs++;
                }
            }
        }

        // Check for fixed lists
        if (nCoarseContribs > 5*fineProcInterface_.interfaceSize())
        {
            WarningIn("processorSAMGInterface::processorSAMGInterface(...)")
                << "Coarse SAMG processor siginificantly bigger than fine: "
                << "nCoarseFaces = " << nCoarseContribs
                << " nFineFaces = " << fineProcInterface_.interfaceSize()
                << endl;
        }

        // Resize arrays to final size
        faceCells_.transfer(dynFaceCells.shrink());
        fineAddressing_.transfer(dynFineAddr.shrink());
        restrictWeights_.transfer(dynRestrictW.shrink());
        restrictAddressing_.transfer(dynRestrictAddr.shrink());

        // Send to slave
        OPstream toNbr(Pstream::blocking, neighbProcNo());

        toNbr<< coarseLevelTable << nCoarseContribs;
    }
    // Slave processor
    else
    {
        // Slave must receive the hash table sent by the master to know how to
        // sort its faceCells_
        IPstream fromNbr(Pstream::blocking, neighbProcNo());

        HashTable<label, labelPair, Hash<labelPair> > masterCoarseLevel =
            HashTable<label, labelPair, Hash<labelPair> >(fromNbr);

        const label nCoarseContribs = readLabel(fromNbr);

        // Filtered prolongation matrix from my side
        crMatrix prolongationT = interfaceProlongation.T();
        const labelList& rRowStart = prolongationT.crAddr().rowStart();
        const labelList& rColumn = prolongationT.crAddr().column();
        const scalarField& rCoeffs = prolongationT.coeffs();
        const label rNRows = prolongationT.crAddr().nRows();

        // Filtered prolongation matrix from neighbour - to obtain restriction,
        // make
        // a transpose
        const labelList& pRowStart =
            nbrInterfaceProlongation.crAddr().rowStart();
        const labelList& pColumn = nbrInterfaceProlongation.crAddr().column();
        const scalarField& pCoeffs = nbrInterfaceProlongation.coeffs();
        const label pNCols = nbrInterfaceProlongation.crAddr().nCols();

        faceCells_.setSize(masterCoarseLevel.size(), -1);
        restrictWeights_.setSize(nCoarseContribs);
        fineAddressing_.setSize(nCoarseContribs);
        restrictAddressing_.setSize(nCoarseContribs);

        labelList coeffMark(pNCols, -1);
        label nCoarseEntries = 0;

        // This loop is only for sorting the faceCells_ on slave side to match
        // the order on master
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
                for
                (
                    label indexP = pRowStart[jr];
                    indexP < pRowStart[jr + 1];
                    indexP++
                )
                {
                    // Column of P, into coarse address
                    label jp = pColumn[indexP];

                    // Array for marking new contributions in the row ir
                    label identify = coeffMark[jp];

                    if (identify == -1)
                    {
                        label address = masterCoarseLevel[labelPair(jp, ir)];

                        // Found a new coarse coeff!
                        identify = address;

                        coeffMark[jp] = address;
                        faceCells_[identify] = ir;
                    }

                    restrictWeights_[nCoarseEntries] =
                        rCoeffs[indexR]*pCoeffs[indexP];
                    fineAddressing_[nCoarseEntries] = jr;
                    restrictAddressing_[nCoarseEntries] = identify;

                    nCoarseEntries++;
                }
            }
        }

    }
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
