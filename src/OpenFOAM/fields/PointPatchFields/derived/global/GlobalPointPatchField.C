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

#include "GlobalPointPatchField.H"
#include "lduMatrix.H"
#include "Map.H"
#include "constraints.H"
#include "PstreamCombineReduceOps.H"

#define OLD_COMBINE_REDUCE 1

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Reduce the field and extract the local values
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
tmp<Field<Type2> > GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::reduceExtractPoint
(
    const tmp<Field<Type2> >& tpField
) const
{
    // Create the global list and insert local values
    if (globalPointPatch_.globalPointSize() > 0)
    {
        // Get addressing
        const labelList& sharedPointAddr =
            globalPointPatch_.sharedPointAddr();

        const Field<Type2>& pField = tpField();

        // Prepare result
        tmp<Field<Type2> > tlpf(new Field<Type2>(sharedPointAddr.size()));
        Field<Type2>& lpf = tlpf();

#       ifdef OLD_COMBINE_REDUCE

        Field<Type2> gpf
        (
            globalPointPatch_.globalPointSize(),
            pTraits<Type2>::zero
        );

        forAll (sharedPointAddr, i)
        {
            gpf[sharedPointAddr[i]] = pField[i];
        }

        combineReduce(gpf, plusEqOp<Field<Type2> >());

        // Extract local data
        forAll (sharedPointAddr, i)
        {
            lpf[i] = gpf[sharedPointAddr[i]];
        }

#       else

        // Pack data into a map
        Map<Type2> dataMap;

        forAll (sharedPointAddr, i)
        {
            dataMap.insert(sharedPointAddr[i], pField[i]);
        }

        // Communicate map
        Pstream::mapCombineGather(dataMap, plusEqOp<Type2>());
        Pstream::mapCombineScatter(dataMap);

        // Extract local data
        forAll (sharedPointAddr, i)
        {
            lpf[i] = dataMap[sharedPointAddr[i]];
        }

#       endif

        return tlpf;
    }
    else
    {
        return tpField;
    }
}


// Reduce the field and extract the local values
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
tmp<Field<Type2> > GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::reduceExtractEdge
(
    const tmp<Field<Type2> >& teField
) const
{
    if (globalPointPatch_.globalEdgeSize() > 0)
    {
        // Bug fix: use map-based communication.  HJ, 18/Nov/2010

        const labelList& sharedEdgeAddr =
            globalPointPatch_.sharedEdgeAddr();

        const Field<Type2>& eField = teField();

        // Prepare result
        tmp<Field<Type2> > tlef(new Field<Type2>(sharedEdgeAddr.size()));
        Field<Type2>& lef = tlef();

#       ifdef OLD_COMBINE_REDUCE

        // Create the global list and insert local values
        Field<Type2> gef
        (
            globalPointPatch_.globalEdgeSize(),
            pTraits<Type2>::zero
        );

        forAll (sharedEdgeAddr, i)
        {
            gef[sharedEdgeAddr[i]] = eField[i];
        }

        combineReduce(gef, plusEqOp<Field<Type2> >());

        // Extract local data
        forAll (sharedEdgeAddr, i)
        {
            lef[i] = gef[sharedEdgeAddr[i]];
        }

#       else

        // Pack data into a map
        Map<Type2> dataMap;

        forAll (sharedEdgeAddr, i)
        {
            dataMap.insert(sharedEdgeAddr[i], eField[i]);
        }

        // Communicate map
        Pstream::mapCombineGather(dataMap, plusEqOp<Type2>());
        Pstream::mapCombineScatter(dataMap);

        // Extract local data
        forAll (sharedEdgeAddr, i)
        {
            lef[i] = dataMap[sharedEdgeAddr[i]];
        }

#       endif

        return tlef;
    }
    else
    {
        return teField;
    }
}


// Add the diagonal/source to the internal field.
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::addFieldTempl
(
    Field<Type2>& pField
) const
{
    // Set the values from the global sum
    tmp<Field<Type2> > trpf =
        reduceExtractPoint<Type2>(patchInternalField(pField));

    Field<Type2>& rpf = trpf();

    // Get addressing
    const labelList& addr = globalPointPatch_.meshPoints();

    forAll (addr, i)
    {
        pField[addr[i]] = rpf[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::GlobalPointPatchField
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
        typename GlobalPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::GlobalPointPatchField
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
        typename GlobalPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(p))
{
    if (!isType<GlobalPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "GlobalPointPatchField<PatchField, Mesh, PointPatch, "
            "GlobalPointPatch, Type>::"
            "GlobalPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index()
            << " not processorPoint type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::GlobalPointPatchField
(
    const GlobalPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        GlobalPointPatch,
        MatrixType,
        Type
    >& ptf,
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
        typename GlobalPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(ptf.patch()))
{
    if (!isType<GlobalPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "GlobalPointPatchField<PatchField, Mesh, PointPatch, "
            "GlobalPointPatch, Type>::"
            "GlobalPointPatchField\n"
            "(\n"
            "    const GlobalPointPatchField<PatchField, Mesh, "
            "PointPatch, GlobalPointPatch, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const PointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::GlobalPointPatchField
(
    const GlobalPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        GlobalPointPatch,
        MatrixType,
        Type
    >& ptf
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename GlobalPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf),
    globalPointPatch_(refCast<const GlobalPointPatch>(ptf.patch()))
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::GlobalPointPatchField
(
    const GlobalPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        GlobalPointPatch,
        MatrixType,
        Type
    >& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename GlobalPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::~GlobalPointPatchField()
{}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (this->isPointField())
    {
        // Create the global list and insert local values
        if (globalPointPatch_.globalPointSize() > 0)
        {
            // Bug fix: use map-based communication.  HJ, 18/Nov/2010

            // Get addressing
            const labelList& sharedPointAddr =
                globalPointPatch_.sharedPointAddr();

            // Get internal field data
            Field<Type> pField = this->patchInternalField();

            // Pack data into a map
            Map<Type> dataMap;

            forAll (sharedPointAddr, i)
            {
                dataMap.insert(sharedPointAddr[i], pField[i]);
            }

            // Communicate map
            // Note: Cannot use reduceExtract, because it uses plusEqOp
            // HJ, 14/Jan/2011
            Pstream::mapCombineGather(dataMap, eqOp<Type>());
            Pstream::mapCombineScatter(dataMap);

            // Extract local data
            Field<Type> lpf(sharedPointAddr.size());

            forAll (sharedPointAddr, i)
            {
                lpf[i] = dataMap[sharedPointAddr[i]];
            }

            // Get addressing and enforce values on all processors
            const labelList& mpAddr = globalPointPatch_.meshPoints();

            // Get internal field to insert values into
            Field<Type>& iF = const_cast<Field<Type>&>(this->internalField());

            forAll (mpAddr, i)
            {
                iF[mpAddr[i]] = lpf[i];
            }
        }
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::addField(Field<Type>& d) const
{
    addFieldTempl(d);
}


// Set boundary condition to matrix
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::setBoundaryCondition
(
    Map<typename MatrixType<Type>::ConstraintType>& fix
) const
{
    // Get addressing
    const labelList& meshPoints = globalPointPatch_.meshPoints();

    forAll (meshPoints, pointI)
    {
        const label curPoint = meshPoints[pointI];

        // Create a constraint. None of the components is fixed
        typename MatrixType<Type>::ConstraintType bc
        (
            curPoint,
            pTraits<Type>::zero,
            pTraits<Type>::zero
        );

        // If not set, add it, otherwise combine with already
        // existing value
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


// Add diagonal coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::addDiag(scalarField& d) const
{
    addFieldTempl(d);
}


// Add source
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::addSource
(
    scalarField& s
) const
{
    addFieldTempl(s);
}


// Add upper/lower coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::addUpperLower
(
    scalarField& eField
) const
{
    // Set the contribution for the local edge coefficients

    // Get the addressing
    const labelList& addr = globalPointPatch_.localEdgeIndices();

    // Get the local elements of the edge field
    tmp<scalarField> tlocalEdgeField(new scalarField(addr.size()));
    scalarField& localEdgeField = tlocalEdgeField();

    forAll (addr, i)
    {
        localEdgeField[i] = eField[addr[i]];
    }

    // Set the edge values
    tmp<scalarField> tref =
        reduceExtractEdge<scalar>(tlocalEdgeField);
    scalarField& ref = tref();

    forAll (addr, i)
    {
        eField[addr[i]] = ref[i];
    }
}


// Get the cut edge coefficients in Amul order
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
tmp<scalarField> GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::cutBouCoeffs
(
    const lduMatrix& m
) const
{
    // Go through all the cut edges.  For all owners pick up the upper
    // and for all the neighbours pick up the lower.

    // Get the indices of cut edges
    const labelList& cutOwn = globalPointPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = globalPointPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = globalPointPatch_.doubleCutEdgeIndices();

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
        cutCoeffs[coeffI] = Upper[cutOwn[edgeI]];
        coeffI++;
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    forAll (cutNei, edgeI)
    {
        cutCoeffs[coeffI] = Lower[cutNei[edgeI]];
        coeffI++;
    }

    // Doubly cut coeficients
    // ~~~~~~~~~~~~~~~~~~~~~~
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
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
tmp<scalarField> GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::cutIntCoeffs
(
    const lduMatrix& m
) const
{
    // Go through all the cut edges.  For all owners pick up the lower
    // and for all the neighbours pick up the upper.

    // Get the indices of cut edges
    const labelList& cutOwn = globalPointPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = globalPointPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = globalPointPatch_.doubleCutEdgeIndices();

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


// Add upper/lower coefficients
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::eliminateUpperLower
(
    scalarField& eField
) const
{
    // Kill the coefficient for cut edges.

    // Get the indices of cut edges
    const labelList& cutOwn = globalPointPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = globalPointPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = globalPointPatch_.doubleCutEdgeIndices();

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


// Complete matrix update on coupled interfaces
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    Mesh,
    PointPatch,
    GlobalPointPatch,
    MatrixType,
    Type
>
::updateInterfaceMatrix
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

    const labelList& mp = globalPointPatch_.meshPoints();

    // Get the multiplication mask for the global side
    const scalarField& cutMask = globalPointPatch_.ownNeiDoubleMask();

    // Get matrix addressing
    const unallocLabelList& L = m.lduAddr().lowerAddr();
    const unallocLabelList& U = m.lduAddr().upperAddr();

    // Note that the addressing is into the local points of the patch.
    // Mesh points is used only for size

    // Coefficients are already ordered in the appropriate way. Just
    // use the counter.
    label coeffI = 0;
    scalarField sumOffDiag(this->size(), 0);

    // Owner side
    // ~~~~~~~~~~
    {
        const labelList& cutOwn = globalPointPatch_.cutEdgeOwnerIndices();
        const labelList& cutOwnStart = globalPointPatch_.cutEdgeOwnerStart();

        forAll (mp, pointI)
        {
            label ownIndex = cutOwnStart[pointI];
            label endOwn = cutOwnStart[pointI + 1];

            for (; ownIndex < endOwn; ownIndex++)
            {
                localMult[pointI] +=
                    cutMask[coeffI]*coeffs[coeffI]
                    *psiInternal[U[cutOwn[ownIndex]]];

                sumOffDiag[pointI] += cutMask[coeffI]*coeffs[coeffI];

                // Multiply the internal side as well
                result[U[cutOwn[ownIndex]]] +=
                    coeffs[coeffI]*psiInternal[mp[pointI]];

                coeffI++;
            }
        }
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    {
        const labelList& cutNei = globalPointPatch_.cutEdgeNeighbourIndices();
        const labelList& cutNeiStart =
            globalPointPatch_.cutEdgeNeighbourStart();

        forAll (mp, pointI)
        {
            label neiIndex = cutNeiStart[pointI];
            label endNei = cutNeiStart[pointI + 1];

            for (; neiIndex < endNei; neiIndex++)
            {
                localMult[pointI] +=
                    cutMask[coeffI]*coeffs[coeffI]
                    *psiInternal[L[cutNei[neiIndex]]];

                sumOffDiag[pointI] += cutMask[coeffI]*coeffs[coeffI];

                // Multiply the internal side as well
                result[L[cutNei[neiIndex]]] +=
                    coeffs[coeffI]*psiInternal[mp[pointI]];

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
    // appropriately.
    {
        const labelList& doubleCut = globalPointPatch_.doubleCutEdgeIndices();

        const labelList& doubleCutOwner = globalPointPatch_.doubleCutOwner();
        const labelList& doubleCutNeighbour =
            globalPointPatch_.doubleCutNeighbour();

        forAll (doubleCut, edgeI)
        {
            // Owner side
            localMult[doubleCutOwner[edgeI]] +=
                cutMask[coeffI]*coeffs[coeffI]*
                psiInternal[U[doubleCut[edgeI]]];

            sumOffDiag[doubleCutOwner[edgeI]] +=
                cutMask[coeffI]*coeffs[coeffI];

            coeffI++;

            // Neighbour side
            localMult[doubleCutNeighbour[edgeI]] +=
                cutMask[coeffI]*coeffs[coeffI]*
                psiInternal[L[doubleCut[edgeI]]];

            sumOffDiag[doubleCutNeighbour[edgeI]] +=
                cutMask[coeffI]*coeffs[coeffI];

            coeffI++;
        }
    }

    // Reduce/extract the result and enforce over all processors

    // Requires global sync points to flush buffers before gather-scatter
    // communications.  Reconsider.  HJ, 29/Mar/2011
    if (Pstream::defaultCommsType == Pstream::nonBlocking)
    {
        IPstream::waitRequests();
        OPstream::waitRequests();
    }

    tmp<Field<scalar> > trpf =
        reduceExtractPoint<scalar>(localMult);

    Field<scalar>& rpf = trpf();

    // Get addressing
    const labelList& addr = globalPointPatch_.meshPoints();

    forAll (addr, i)
    {
        result[addr[i]] += rpf[i];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
