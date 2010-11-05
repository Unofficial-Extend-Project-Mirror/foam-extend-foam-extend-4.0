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

#include "processorFvPatchVectorNFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#define VectorNMatrixInterfaceFunc(Type)                                        \
template <>                                                                     \
void processorFvPatchField<Type>::initInterfaceMatrixUpdate                     \
(                                                                               \
    const Field<Type>& psiInternal,                                             \
    Field<Type>&,                                                               \
    const BlockLduMatrix<Type>&,                                                \
    const CoeffField<Type>&,                                                    \
    const Pstream::commsTypes commsType                                         \
) const                                                                         \
{                                                                               \
    procPatch_.compressedSend                                                   \
    (                                                                           \
        commsType,                                                              \
        this->patch().patchInternalField(psiInternal)()                         \
    );                                                                          \
}                                                                               \
                                                                                \
template <>                                                                     \
void processorFvPatchField<Type>::updateInterfaceMatrix                         \
(                                                                               \
    const Field<Type>& psiInternal,                                             \
    Field<Type>& result,                                                        \
    const BlockLduMatrix<Type>&,                                                \
    const CoeffField<Type>& coeffs,                                             \
    const Pstream::commsTypes commsType                                         \
) const                                                                         \
{                                                                               \
    Field<Type> pnf(this->size());                                              \
                                                                                \
    if (coeffs.activeType() == blockCoeffBase::SCALAR)                          \
    {                                                                           \
        pnf = coeffs.asScalar() *                                               \
            procPatch_.compressedReceive<Type>(commsType, this->size())();      \
    }                                                                           \
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)                     \
    {                                                                           \
        pnf = cmptMultiply(coeffs.asLinear(),                                   \
            procPatch_.compressedReceive<Type>(commsType, this->size())()       \
        );                                                                      \
    }                                                                           \
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)                     \
    {                                                                           \
        pnf = coeffs.asSquare() &                                               \
            procPatch_.compressedReceive<Type>(commsType, this->size())();      \
    }                                                                           \
                                                                                \
    const unallocLabelList& faceCells = this->patch().faceCells();              \
                                                                                \
    forAll(faceCells, facei)                                                    \
    {                                                                           \
        result[faceCells[facei]] -= pnf[facei];                                 \
    }                                                                           \
}


#define doMakePatchTypeField(type, Type, args...)                               \
    VectorNMatrixInterfaceFunc(type)                                            \
                                                                                \
    makePatchTypeField(fvPatch##Type##Field, processorFvPatch##Type##Field);

forAllVectorNTypes(doMakePatchTypeField)

#undef doMakePatchTypeField

#undef VectorNMatrixInterfaceFunc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
