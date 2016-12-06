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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "transformField.H"
#include "coeffFields.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const processorFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(0),
    receiveBuf_(0),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (debug && !ptf.ready())
    {
        FatalErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const processorFvPatchField<Type>& ptf,\n"
            "    const DimensionedField<Type, volMesh>& iF\n"
            ")"
        )   << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::~processorFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::processorFvPatchField<Type>::patchNeighbourField() const
{
    if (debug && !this->ready())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> >"
            "processorFvPatchField<Type>::patchNeighbourField() const"
        )   << "On patch " << procPatch_.name()
            << " outstanding request."
            << abort(FatalError);
    }

    // Warning: returning own patch field, which only after update stores
    // actual neighbour data
    // HJ, 14/May/2009
    return *this;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        // Collect data into send buffer
        sendBuf_ = this->patchInternalField();

        if (commsType == Pstream::nonBlocking)
        {
            // Fast path. Receive into *this

            outstandingRecvRequest_ = Pstream::nRequests();

            IPstream::read
            (
                Pstream::nonBlocking,
                procPatch_.neighbProcNo(),
                reinterpret_cast<char*>(this->begin()),
                this->byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );

            outstandingSendRequest_ = Pstream::nRequests();

            OPstream::write
            (
                Pstream::nonBlocking,
                procPatch_.neighbProcNo(),
                reinterpret_cast<const char*>(sendBuf_.begin()),
                this->byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }
        else
        {
            procPatch_.send(commsType, sendBuf_);
        }
    }
}


template<class Type>
void Foam::processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if (commsType == Pstream::nonBlocking)
        {
            // Fast path. Received into *this

            if
            (
                outstandingRecvRequest_ >= 0
             && outstandingRecvRequest_ < Pstream::nRequests()
            )
            {
                Pstream::waitRequest(outstandingRecvRequest_);
            }
            outstandingSendRequest_ = -1;
            outstandingRecvRequest_ = -1;
        }
        else
        {
            procPatch_.receive<Type>(commsType, *this);
        }

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), *this);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::processorFvPatchField<Type>::snGrad() const
{
    return this->patch().deltaCoeffs()*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
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
    scalarSendBuf_ = this->patch().patchInternalField(psiInternal);

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorIn
            (
                "void processorFvPatchField<Type>::initInterfaceMatrixUpdate\n"
                "(\n"
                "    const scalarField& psiInternal,\n"
                "    scalarField&,\n"
                "    const lduMatrix&,\n"
                "    const scalarField&,\n"
                "    const direction,\n"
                "    const Pstream::commsTypes commsType,\n"
                "    const bool switchToLhs\n"
                ") const"
            )   << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        // Check buffer size
        scalarReceiveBuf_.setSize(this->size());

        outstandingRecvRequest_ = Pstream::nRequests();

        IPstream::read
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = Pstream::nRequests();

        OPstream::write
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.send(commsType, scalarSendBuf_);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
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
        scalarReceiveBuf_.setSize(this->size());

        procPatch_.receive<scalar>(commsType, scalarReceiveBuf_);
    }

    // The data is now in scalarReceiveBuf_ for both cases

    // Transform according to the transformation tensor
    transformCoupleField(scalarReceiveBuf_, cmpt);

    // Multiply the field by coefficients and add into the result

    const unallocLabelList& faceCells = this->patch().faceCells();

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

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const Field<Type>& psiInternal,
    Field<Type>&,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    sendBuf_ = this->patch().patchInternalField(psiInternal);

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorIn
            (
                "void processorFvPatchField<Type>::initInterfaceMatrixUpdate\n"
                "(\n"
                "    const Field<Type>& psiInternal,\n"
                "    Field<Type>&,\n"
                "    const BlockLduMatrix<Type>&,\n"
                "    const CoeffField<Type>&,\n"
                "    const Pstream::commsTypes commsType,\n"
                "    const bool switchToLhs\n"
                ") const"
            )   << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        // Check buffer size
        receiveBuf_.setSize(this->size());

        outstandingRecvRequest_ = Pstream::nRequests();

        IPstream::read
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(receiveBuf_.begin()),
            receiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = Pstream::nRequests();

        OPstream::write
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(sendBuf_.begin()),
            sendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.send
        (
            commsType,
            this->patch().patchInternalField(psiInternal)()
        );
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
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
        receiveBuf_.setSize(this->size());

        procPatch_.receive<Type>(commsType, receiveBuf_);
    }

    // The data is now in receiveBuf_ for both cases

    // Multiply neighbour field with coeffs and re-use buffer for result
    // of multiplication
    multiply(receiveBuf_, coeffs, receiveBuf_);

    const unallocLabelList& faceCells = this->patch().faceCells();

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

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
bool Foam::processorFvPatchField<Type>::ready() const
{
    if
    (
        outstandingSendRequest_ >= 0
     && outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished = Pstream::finishedRequest(outstandingSendRequest_);

        if (!finished)
        {
            return false;
        }
    }
    outstandingSendRequest_ = -1;

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished = Pstream::finishedRequest(outstandingRecvRequest_);

        if (!finished)
        {
            return false;
        }
    }
    outstandingRecvRequest_ = -1;

    return true;
}


// ************************************************************************* //
