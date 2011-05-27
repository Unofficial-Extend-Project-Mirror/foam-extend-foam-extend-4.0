/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "ProcessorPointPatchField.H"
#include "lduMatrix.H"
#include "Map.H"
#include "transformField.H"
#include "constraints.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::resizeBuf
(
    List<char>& buf, 
    const label size
) const
{
    if (buf.size() < size)
    {
        buf.setSize(size);
    }
}


// Raw field sending and receiving
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
sendField
(
    const tmp<Field<Type2> >& tf,
    const Pstream::commsTypes commsType
) const
{
    const Field<Type2>& f = tf();

    //HJ: This needs complete rewrite:
    // - move communications into a patch
    // - allow for various types of communication
    // HJ, 15/Apr/2009

    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        OPstream::write
        (
            commsType,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(f.begin()),
            f.byteSize()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        resizeBuf(receiveBuf_, f.size()*sizeof(Type));

        IPstream::read
        (
            commsType,
            procPatch_.neighbProcNo(),
            receiveBuf_.begin(),
            receiveBuf_.size()
        );

        resizeBuf(sendBuf_, f.byteSize());
        memcpy(sendBuf_.begin(), f.begin(), f.byteSize());

        OPstream::write
        (
            commsType,
            procPatch_.neighbProcNo(),
            sendBuf_.begin(),
            f.byteSize()
        );
    }
    else
    {
        FatalErrorIn("ProcessorPointPatchField::send")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }

    // Not using non-blocking comms
//     if (commsType == Pstream::nonBlocking)
//     {
//         FatalErrorIn("void ProcessorPointPatchField::sendField")
//             << "Non-blocking comms not implemented"
//             << abort(FatalError);
//     }

//     OPstream::write
//     (
//         commsType,
//         procPatch_.neighbProcNo(),
//         reinterpret_cast<const char*>(f.begin()),
//         f.byteSize()
//     );

    tf.clear();
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
tmp<Field<Type2> >
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
receivePointField
(
    const Pstream::commsTypes commsType
) const
{
    tmp<Field<Type2> > tf(new Field<Type2>(this->size()));

    IPstream::read
    (
        commsType,
        procPatch_.neighbProcNo(),
        reinterpret_cast<char*>(tf().begin()),
        tf().byteSize()
    );

    return tf;
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
tmp<Field<Type2> >
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
receiveEdgeField
(
    const Pstream::commsTypes commsType
) const
{
    tmp<Field<Type2> > tf
    (
        new Field<Type2>(procPatch_.localEdgeIndices().size())
    );

    IPstream::read
    (
        commsType,
        procPatch_.neighbProcNo(),
        reinterpret_cast<char*>(tf().begin()),
        tf().byteSize()
    );

    return tf;
}


// Initialise diagonal/source update.
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initAddFieldTempl
(
    const Pstream::commsTypes commsType,
    const Field<Type2>& pField
) const
{
    sendField(patchInternalField(pField), commsType);
}


// Add the diagonal/source to the internal field.
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
addFieldTempl
(
    const Pstream::commsTypes commsType,
    Field<Type2>& pField
) const
{
    // Get the neighbour side values
    tmp<Field<Type2> > tpNeighbour = receivePointField<Type2>(commsType);
    addToInternalField(pField, tpNeighbour());
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
ProcessorPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p)),
    sendBuf_(),
    receiveBuf_()
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
ProcessorPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p)),
    sendBuf_(),
    receiveBuf_()
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
    <PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>& ptf,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch())),
    sendBuf_(),
    receiveBuf_()
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
    <PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>& ptf
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch())),
    sendBuf_(),
    receiveBuf_()
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
    <PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch())),
    sendBuf_(),
    receiveBuf_()
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
~ProcessorPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Initialise field transfer
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if (this->isPointField())
        {
            initAddFieldTempl(Pstream::blocking, this->internalField());
        }
    }
}


// Initialise field transfer
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if (this->isPointField())
        {
            // Get the neighbour side values
            tmp<Field<Type> > tpNeighbour = receivePointField<Type>(commsType);
            Field<Type>& tpn = tpNeighbour();

            if (doTransform())
            {
                const processorPolyPatch& ppp = procPatch_.procPolyPatch();
                const tensorField& forwardT = ppp.forwardT();

                transform(tpn, forwardT[0], tpn);
            }

            // Average over two sides
            tpn = 0.5*(patchInternalField(this->internalField()) + tpn);

            // Get internal field to insert values into
            Field<Type>& iF = const_cast<Field<Type>&>(this->internalField());

            this->setInInternalField(iF, tpn);
        }
    }
}


// Initialise field transfer
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initAddField() const
{
    initAddFieldTempl(Pstream::blocking, this->internalField());
}


// Add field
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
addField(Field<Type>& f) const
{
    addFieldTempl(Pstream::blocking, f);
}


// Set boundary condition to matrix
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
setBoundaryCondition
(
    Map<typename MatrixType<Type>::ConstraintType> & fix
) const
{
    // get addressing
    const labelList& meshPoints = procPatch_.meshPoints();

    forAll (meshPoints, pointI)
    {
        const label curPoint = meshPoints[pointI];

        // create a constraint.  None of the components is fixed
        typename MatrixType<Type>::ConstraintType bc
        (
            curPoint,
            pTraits<Type>::zero,
            pTraits<Type>::zero
        );

        // If pointer is not set, add it, otherwise combine with
        // already existing value
        if (!fix.found(curPoint))
        {
            fix.insert(curPoint, bc);
        }
        else
        {
            fix[curPoint].combine(bc);
        }
    }
}


// Initialise transfer of diagonal coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initAddDiag(const scalarField& d) const
{
    initAddFieldTempl(Pstream::blocking, d);
}


// Initialise transfer of source
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initAddSource(const scalarField& s) const
{
    initAddFieldTempl(Pstream::blocking, s);
}


// Add diagonal coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
addDiag(scalarField& d) const
{
    addFieldTempl(Pstream::blocking, d);
}


// Add source
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
addSource(scalarField& s) const
{
    addFieldTempl(Pstream::blocking, s);
}


// Initialise transfer of upper/lower coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initAddUpperLower
(
    const scalarField& eField
) const
{
    // Gather the data from the given field and send it. Note that the
    // size of the field is equal to the number of edges on the field
    // and NOT the number of points.  HJ, 14/Nov/2001

    // Get the addressing
    const labelList& me = procPatch_.localEdgeIndices();

    tmp<scalarField> tresult(new scalarField(me.size()));
    scalarField& result = tresult();

    forAll (me, edgeI)
    {
        result[edgeI] = eField[me[edgeI]];
    }

    // Send the result
    sendField(tresult);
}


// Add upper/lower coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
addUpperLower
(
    scalarField& eField
) const
{
    // Add the contribution for the local edge coefficients

    // Get the neighbour side values
    tmp<scalarField> teNeighbour = receiveEdgeField<scalar>();
    scalarField& eNeighbour = teNeighbour();

    // Get the addressing
    const labelList& me = procPatch_.localEdgeIndices();

    forAll (me, edgeI)
    {
        eField[me[edgeI]] += eNeighbour[edgeI];
    }
}


// Get the cut edge coefficients in Amul order
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
tmp<scalarField>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
cutBouCoeffs
(
    const lduMatrix& m
) const
{
    // Go through all the cut edges.  For all owners pick up the upper
    // and for all the neighbours pick up the lower.

    // Get the indices of cut edges
    const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

    // Get matrix coefficients
    const scalarField& Lower = m.lower();
    const scalarField& Upper = m.upper();

    tmp<scalarField> tcutCoeffs
    (
        new scalarField(cutOwn.size() + cutNei.size() + 2*doubleCut.size(), 0)
    );
    scalarField& cutCoeffs = tcutCoeffs();

    label coeffI = 0;
        
    // Owner side

    forAll (cutOwn, edgeI)
    {
        cutCoeffs[coeffI] = Upper[cutOwn[edgeI]];
        coeffI++;
    }

    // Neighbour side

    forAll (cutNei, edgeI)
    {
        cutCoeffs[coeffI] = Lower[cutNei[edgeI]];
        coeffI++;
    }

    // Doubly cut coeficients

    forAll (doubleCut, edgeI)
    {
        cutCoeffs[coeffI] = Upper[doubleCut[edgeI]];
        coeffI++;

        cutCoeffs[coeffI] = Lower[doubleCut[edgeI]];
        coeffI++;
    }

    return tcutCoeffs;
}


// Get the cut edge coefficients in Tmul order
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
tmp<scalarField>
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
cutIntCoeffs
(
    const lduMatrix& m
) const
{
    // Go through all the cut edges.  For all owners pick up the lower
    // and for all the neighbours pick up the upper.

    // Get the indices of cut edges
    const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

    // Get matrix coefficients
    const scalarField& Lower = m.lower();
    const scalarField& Upper = m.upper();

    tmp<scalarField> tcutCoeffs
    (
        new scalarField(cutOwn.size() + cutNei.size() + 2*doubleCut.size(), 0)
    );
    scalarField& cutCoeffs = tcutCoeffs();

    label coeffI = 0;

    // Owner side
    // ~~~~~~~~~~
    forAll (cutOwn, edgeI)
    {
        cutCoeffs[coeffI] = Lower[cutOwn[edgeI]];
        coeffI++;
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    forAll (cutNei, edgeI)
    {
        cutCoeffs[coeffI] = Upper[cutNei[edgeI]];
        coeffI++;
    }

    // Doubly cut coeficients
    // ~~~~~~~~~~~~~~~~~~~~~~
    forAll (doubleCut, edgeI)
    {
        cutCoeffs[coeffI] = Lower[doubleCut[edgeI]];
        coeffI++;

        cutCoeffs[coeffI] = Upper[doubleCut[edgeI]];
        coeffI++;
    }

    return tcutCoeffs;
}


// Eliminate upper/lower coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
eliminateUpperLower
(
    scalarField& eField
) const
{
    // Kill the coefficient for cut edges.

    // Get the indices of cut edges
    const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

    // Owner side
    // ~~~~~~~~~~
    forAll (cutOwn, edgeI)
    {
        eField[cutOwn[edgeI]] = 0;
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    forAll (cutNei, edgeI)
    {
        eField[cutNei[edgeI]] = 0;
    }

    // Doubly cut edges
    // ~~~~~~~~~~~~~~~~
    forAll (doubleCut, edgeI)
    {
        eField[doubleCut[edgeI]] = 0;
    }
}


// Initialise matrix update on coupled interfaces
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    tmp<scalarField> tlocalMult(new scalarField(this->size(), 0));
    scalarField& localMult = tlocalMult();

    const labelList& mp = procPatch_.meshPoints();

    // Get matrix addressing
    const unallocLabelList& L = m.lduAddr().lowerAddr();
    const unallocLabelList& U = m.lduAddr().upperAddr();

    // Note that the addressing is into the local points of the patch.
    // Mesh points is used only for size

    // Get the multiplication mask to exclude all unwanted local multiplies.
    // An example of this is an internal edge between two points which
    // belong to two different processor patches
    const scalarField& cutMask = procPatch_.ownNeiDoubleMask();

    // Coefficients are already ordered in the appropriate way. Just
    // use the counter.  
    label coeffI = 0;

    // Owner side
    // ~~~~~~~~~~
    {
        const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
        const labelList& cutOwnStart = procPatch_.cutEdgeOwnerStart();

        forAll (mp, pointI)
        {
            label ownIndex = cutOwnStart[pointI];
            label endOwn = cutOwnStart[pointI + 1];

            for (; ownIndex < endOwn; ownIndex++)
            {
                localMult[pointI] +=
                    coeffs[coeffI]*psiInternal[U[cutOwn[ownIndex]]];

                // Multiply the internal side as well using the cut mask
                result[U[cutOwn[ownIndex]]] +=
                    cutMask[coeffI]*coeffs[coeffI]*psiInternal[mp[pointI]];

                coeffI++;
            }
        }
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    {
        const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
        const labelList& cutNeiStart = procPatch_.cutEdgeNeighbourStart();

        forAll (mp, pointI)
        {
            label neiIndex = cutNeiStart[pointI];
            label endNei = cutNeiStart[pointI + 1];

            for (; neiIndex < endNei; neiIndex++)
            {
                localMult[pointI] +=
                    coeffs[coeffI]*psiInternal[L[cutNei[neiIndex]]];

                // Multiply the internal side as well using the cut mask
                result[L[cutNei[neiIndex]]] +=
                    cutMask[coeffI]*coeffs[coeffI]*psiInternal[mp[pointI]];

                coeffI++;
            }
        }
    }

    // Doubly cut coefficients
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // There exists a possibility of having an internal edge for a
    // point on the processor patch which is in fact connected to
    // another point of the same patch.  This particular nastiness
    // introduces a deformation in the solution because the edge is
    // either multiplied twice or not at all.  For this purpose, the
    // offending edges need to be separated out and multiplied
    // appropriately.  This will only happen for cell tetrahedral
    // decomposition and is generally nasty.  
    // No need for cut mask here
    {
        const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

        const labelList& doubleCutOwner = procPatch_.doubleCutOwner();
        const labelList& doubleCutNeighbour = procPatch_.doubleCutNeighbour();

        forAll (doubleCut, edgeI)
        {
            // Owner side
            localMult[doubleCutOwner[edgeI]] +=
                coeffs[coeffI]*psiInternal[U[doubleCut[edgeI]]];
            coeffI++;

            // Neighbour side
            localMult[doubleCutNeighbour[edgeI]] +=
                coeffs[coeffI]*psiInternal[L[doubleCut[edgeI]]];
            coeffI++;
        }
    }

    // Add the local multiplication to this side as well

    forAll (mp, pointI)
    {
        result[mp[pointI]] += localMult[pointI];
    }

    // Send the localMult
    sendField(tlocalMult, commsType);
}


// Complete matrix update on coupled interfaces
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
void
ProcessorPointPatchField
<PatchField, Mesh, PointPatch, ProcessorPointPatch, MatrixType, Type>::
updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    // Get the neighbour side multiplication
    tmp<scalarField> tneiMult = receivePointField<scalar>(commsType);
    this->addToInternalField(result, tneiMult());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
