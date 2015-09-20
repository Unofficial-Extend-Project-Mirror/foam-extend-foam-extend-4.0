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

#include "cyclicFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{
    if (!isType<cyclicFaPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicFaPatchField<Type>::cyclicFaPatchField\n"
            "(\n"
            "    const cyclicFaPatchField<Type>& ptf,\n"
            "    const faPatch& p,\n"
            "    const DimensionedField<Type, areaMesh>& iF,\n"
            "    const faPatchFieldMapper& mapper\n"
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
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    coupledFaPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{
    if (!isType<cyclicFaPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicFaPatchField<Type>::cyclicFaPatchField\n"
            "(\n"
            "    const faPatch& p,\n"
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

    this->evaluate(Pstream::blocking);
}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > cyclicFaPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();
    const unallocLabelList& faceCells = cyclicPatch_.faceCells();

    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();

    label sizeby2 = this->size()/2;

    if (doTransform())
    {
        for (label facei=0; facei<sizeby2; facei++)
        {
            pnf[facei] = transform
            (
                forwardT()[0], iField[faceCells[facei + sizeby2]]
            );

            pnf[facei + sizeby2] = transform
            (
                reverseT()[0], iField[faceCells[facei]]
            );
        }
    }
    else
    {
        for (label facei=0; facei<sizeby2; facei++)
        {
            pnf[facei] = iField[faceCells[facei + sizeby2]];
            pnf[facei + sizeby2] = iField[faceCells[facei]];
        }
    }

    return tpnf;
}


template<class Type>
void cyclicFaPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    scalarField pnf(this->size());

    label sizeby2 = this->size()/2;
    const unallocLabelList& faceCells = cyclicPatch_.faceCells();

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
