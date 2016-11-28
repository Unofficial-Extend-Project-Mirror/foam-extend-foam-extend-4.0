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

#include "ProcessorBlockAMGInterfaceField.H"
#include "processorLduInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::ProcessorBlockAMGInterfaceField<Type>::ProcessorBlockAMGInterfaceField
(
    const AMGInterface& AMGCp,
    const BlockLduInterfaceField<Type>& fineInterfaceField
)
:
    BlockAMGInterfaceField<Type>(AMGCp, fineInterfaceField),
    procInterface_(refCast<const processorAMGInterface>(AMGCp)),
    doTransform_(false),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0)
{
    // If the interface based on a patch this must be taken care specially of
    if (isA<ProcessorBlockLduInterfaceField<Type> >(fineInterfaceField))
    {
        const ProcessorBlockLduInterfaceField<Type>& p =
            refCast<const ProcessorBlockLduInterfaceField<Type> >
            (
                fineInterfaceField
            );

        doTransform_ = p.doTransform();
    }
    else if (isA<processorLduInterfaceField>(fineInterfaceField))
    {
        const processorLduInterfaceField& p =
            refCast<const processorLduInterfaceField >(fineInterfaceField);

        doTransform_ = p.doTransform();
    }
    else
    {
        FatalErrorIn("ProcessorBlockAMGInterfaceField<Type> Constructor")
            << "fineInterface must be of processor type and either" << endl
            << "    ProcessorBlockLduInterfaceField<Type> or " << endl
            << "    processorFvPatchField<Type> " << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::ProcessorBlockAMGInterfaceField<Type>::
~ProcessorBlockAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::ProcessorBlockAMGInterfaceField<Type>::initInterfaceMatrixUpdate
(
    const Field<Type>& psiInternal,
    Field<Type>&,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    label oldWarn = Pstream::warnComm;
    Pstream::warnComm = comm();

    sendBuf_ = procInterface_.interfaceInternalField(psiInternal);

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        receiveBuf_.setSize(sendBuf_.size());
        outstandingRecvRequest_ = Pstream::nRequests();
        IPstream::read
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<char*>(receiveBuf_.begin()),
            receiveBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );

        outstandingSendRequest_ = Pstream::nRequests();
        OPstream::write
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<const char*>(sendBuf_.begin()),
            sendBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );
    }
    else
    {
        procInterface_.send
        (
            commsType,
            procInterface_.interfaceInternalField(psiInternal)()
        );
    }

    // Mark as ready for update
    const_cast<ProcessorBlockAMGInterfaceField<Type>&>(*this).updatedMatrix() =
        false;

    Pstream::warnComm = oldWarn;
}


template<class Type>
void Foam::ProcessorBlockAMGInterfaceField<Type>::updateInterfaceMatrix
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            Pstream::waitRequest(outstandingRecvRequest_);
        }

        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;
    }
    else
    {
        // Check size
        receiveBuf_.setSize(sendBuf_.size());

        procInterface_.receive<Type>(commsType, receiveBuf_);
    }

    // The data is now in receiveBuf_ for both cases

    // Transformation missing.  Is it needed?  HJ, 28/Nov/2016

    // Multiply neighbour field with coeffs and re-use buffer for result
    // of multiplication
    multiply(receiveBuf_, coeffs, receiveBuf_);

    const unallocLabelList& faceCells = procInterface_.faceCells();

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += receiveBuf_[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= receiveBuf_[elemI];
        }
    }

    // Mark as updated
    const_cast<ProcessorBlockAMGInterfaceField<Type>&>(*this).updatedMatrix() =
        true;
}


// ************************************************************************* //
