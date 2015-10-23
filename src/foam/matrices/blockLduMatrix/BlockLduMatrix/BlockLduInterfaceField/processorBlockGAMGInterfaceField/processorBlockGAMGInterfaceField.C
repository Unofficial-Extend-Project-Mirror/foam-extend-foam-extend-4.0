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

#include "processorBlockGAMGInterfaceField.H"
#include "processorLduInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorBlockGAMGInterfaceField<Type>::processorBlockGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const BlockLduInterfaceField<Type>& fineInterfaceField
)
:
    BlockGAMGInterfaceField<Type>(GAMGCp, fineInterfaceField),
    procInterface_(refCast<const processorGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    // If the interface based on a patch this must be taken care specially of
    if (isA<processorBlockLduInterfaceField<Type> >(fineInterfaceField))
    {
        const processorBlockLduInterfaceField<Type>& p =
            refCast<const processorBlockLduInterfaceField<Type> >
            (
                fineInterfaceField
            );

        doTransform_ = p.doTransform();
        rank_ = p.rank();
    }
    else if (isA<processorLduInterfaceField>(fineInterfaceField))
    {
        const processorLduInterfaceField& p =
            refCast<const processorLduInterfaceField >(fineInterfaceField);

        doTransform_ = p.doTransform();
        rank_ = p.rank();
    }
    else
    {
        FatalErrorIn("processorBlockGAMGInterfaceField<Type> Constructor")
            << "fineInterface must be of processor type and either" << endl
            << "    processorBlockLduInterfaceField<Type> or " << endl
            << "    processorFvPatchField<Type> " << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::processorBlockGAMGInterfaceField<Type>::
~processorBlockGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorBlockGAMGInterfaceField<Type>::initInterfaceMatrixUpdate
(
    const Field<Type>& psiInternal,
    Field<Type>&,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    procInterface_.compressedSend
    (
        commsType,
        procInterface_.interfaceInternalField(psiInternal)()
    );
}

template<class Type>
void Foam::processorBlockGAMGInterfaceField<Type>::updateInterfaceMatrix
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    Field<Type> pnf
    (
        coeffs.size()
    );

    if (coeffs.activeType() == blockCoeffBase::SCALAR)
    {
        pnf = coeffs.asScalar()*
            procInterface_.compressedReceive<Type>(commsType, this->size())();
    }
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)
    {
        pnf = cmptMultiply
        (
            coeffs.asLinear(),
            procInterface_.compressedReceive<Type>(commsType, this->size())()
        );
    }
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)
    {
        pnf = coeffs.asSquare() &
            procInterface_.compressedReceive<Type>(commsType, this->size())();
    }

    const unallocLabelList& faceCells = procInterface_.faceCells();

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= pnf[elemI];
        }
    }
}


// ************************************************************************* //
