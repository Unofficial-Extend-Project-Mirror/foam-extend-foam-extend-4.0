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

#include "processorAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        AMGInterfaceField,
        processorAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorAMGInterfaceField::processorAMGInterfaceField
(
    const AMGInterface& AMGCp,
    const lduInterfaceField& fineInterfaceField
)
:
    AMGInterfaceField(AMGCp, fineInterfaceField),
    procInterface_(refCast<const processorAMGInterface>(AMGCp)),
    doTransform_(false),
    rank_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    const processorLduInterfaceField& p =
        refCast<const processorLduInterfaceField>(fineInterfaceField);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorAMGInterfaceField::~processorAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorAMGInterfaceField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    label oldWarn = Pstream::warnComm;
    Pstream::warnComm = comm();

    scalarSendBuf_ = procInterface_.interfaceInternalField(psiInternal);

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        scalarReceiveBuf_.setSize(scalarSendBuf_.size());
        outstandingRecvRequest_ = Pstream::nRequests();
        IPstream::read
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );

        outstandingSendRequest_ = Pstream::nRequests();
        OPstream::write
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
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
    const_cast<processorAMGInterfaceField&>(*this).updatedMatrix() = false;

    Pstream::warnComm = oldWarn;
}


void Foam::processorAMGInterfaceField::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
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
        scalarReceiveBuf_.setSize(scalarSendBuf_.size());

        procInterface_.receive<scalar>(commsType, scalarReceiveBuf_);
    }

    // The data is now in scalarReceiveBuf_ for both cases

    // Transform according to the transformation tensor
    transformCoupleField(scalarReceiveBuf_, cmpt);

    // Multiply the field by coefficients and add into the result

    const unallocLabelList& faceCells = procInterface_.faceCells();

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*scalarReceiveBuf_[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*scalarReceiveBuf_[elemI];
        }
    }

    // Mark as updated
    const_cast<processorAMGInterfaceField&>(*this).updatedMatrix() = true;
}


// ************************************************************************* //
