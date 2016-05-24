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

#include "ProcessorBlockAMGInterfaceField.H"
#include "processorLduInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


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
    doTransform_(false)
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
    procInterface_.compressedSend
    (
        commsType,
        procInterface_.interfaceInternalField(psiInternal)()
    );
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
    Field<Type> pnf
    (
        procInterface_.compressedReceive<Type>(commsType, this->size())
    );

    // Multiply neighbour field with coeffs and re-use pnf for result
    // of multiplication
    multiply(pnf, coeffs, pnf);

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
