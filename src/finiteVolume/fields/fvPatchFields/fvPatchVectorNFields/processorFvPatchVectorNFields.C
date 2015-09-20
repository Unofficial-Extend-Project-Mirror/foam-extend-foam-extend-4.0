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

#include "processorFvPatchVectorNFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#define VectorNMatrixInterfaceFunc(Type)                                      \
template <>                                                                   \
void processorFvPatchField<Type>::initInterfaceMatrixUpdate                   \
(                                                                             \
    const Field<Type>& psiInternal,                                           \
    Field<Type>&,                                                             \
    const BlockLduMatrix<Type>&,                                              \
    const CoeffField<Type>&,                                                  \
    const Pstream::commsTypes commsType,                                      \
    const bool switchToLhs                                                    \
) const                                                                       \
{                                                                             \
    procPatch_.compressedSend                                                 \
    (                                                                         \
        commsType,                                                            \
        this->patch().patchInternalField(psiInternal)()                       \
    );                                                                        \
}                                                                             \
                                                                              \
template <>                                                                   \
void processorFvPatchField<Type>::updateInterfaceMatrix                       \
(                                                                             \
    const Field<Type>& psiInternal,                                           \
    Field<Type>& result,                                                      \
    const BlockLduMatrix<Type>&,                                              \
    const CoeffField<Type>& coeffs,                                           \
    const Pstream::commsTypes commsType,                                      \
    const bool switchToLhs                                                    \
) const                                                                       \
{                                                                             \
    Field<Type> pnf(this->size());                                            \
                                                                              \
    if (coeffs.activeType() == blockCoeffBase::SCALAR)                        \
    {                                                                         \
        pnf = coeffs.asScalar() *                                             \
            procPatch_.compressedReceive<Type>(commsType, this->size())();    \
    }                                                                         \
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)                   \
    {                                                                         \
        pnf = cmptMultiply(coeffs.asLinear(),                                 \
            procPatch_.compressedReceive<Type>(commsType, this->size())()     \
        );                                                                    \
    }                                                                         \
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)                   \
    {                                                                         \
        pnf = coeffs.asSquare() &                                             \
            procPatch_.compressedReceive<Type>(commsType, this->size())();    \
    }                                                                         \
                                                                              \
    const unallocLabelList& faceCells = this->patch().faceCells();            \
                                                                              \
    if (switchToLhs)                                                          \
    {                                                                         \
        forAll(faceCells, elemI)                                              \
        {                                                                     \
            result[faceCells[elemI]] += pnf[elemI];                           \
        }                                                                     \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        forAll(faceCells, elemI)                                              \
        {                                                                     \
            result[faceCells[elemI]] -= pnf[elemI];                           \
        }                                                                     \
    }                                                                         \
}


#define doMakePatchTypeField(type, Type, args...)                             \
                                                                              \
VectorNMatrixInterfaceFunc(type)                                              \
                                                                              \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatch##Type##Field,                                                     \
    processorFvPatch##Type##Field                                             \
);

forAllVectorNTypes(doMakePatchTypeField)

#undef doMakePatchTypeField

#undef VectorNMatrixInterfaceFunc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
