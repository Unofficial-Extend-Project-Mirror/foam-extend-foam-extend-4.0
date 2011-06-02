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

#include "cyclicFvPatchVectorNFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#define VectorNMatrixInterfaceFunc(Type)                                        \
template<>																		\
void cyclicFvPatchField<Type>::updateInterfaceMatrix							\
(																				\
    const Field<Type>& psiInternal,                                             \
    Field<Type>& result,                                                        \
    const BlockLduMatrix<Type>&,                                                \
    const CoeffField<Type>& coeffs,                                             \
    const Pstream::commsTypes commsType                                         \
) const																			\
{																				\
    Field<Type> pnf(this->size());												\
																				\
    label sizeby2 = this->size()/2;												\
    const unallocLabelList& faceCells = cyclicPatch_.faceCells();				\
																				\
    for (label facei=0; facei<sizeby2; facei++)									\
    {																			\
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];					\
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];					\
    }																			\
																				\
    if (coeffs.activeType() == blockCoeffBase::SCALAR)                          \
    {                                                                           \
        pnf = coeffs.asScalar() * pnf;                                          \
    }                                                                           \
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)                     \
    {                                                                           \
        pnf = cmptMultiply(coeffs.asLinear(), pnf);								\
    }                                                                           \
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)                     \
    {                                                                           \
        pnf = coeffs.asSquare() & pnf;											\
    }                                                                           \
                                                                                \
    forAll(faceCells, elemI)													\
    {																			\
        result[faceCells[elemI]] -= pnf[elemI];									\
    }																			\
}


#define doMakePatchTypeField(type, Type, args...)                               \
    VectorNMatrixInterfaceFunc(type)                                            \
                                                                                \
    makePatchTypeField(fvPatch##Type##Field, cyclicFvPatch##Type##Field);

forAllVectorNTypes(doMakePatchTypeField)

#undef doMakePatchTypeField

#undef VectorNMatrixInterfaceFunc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
