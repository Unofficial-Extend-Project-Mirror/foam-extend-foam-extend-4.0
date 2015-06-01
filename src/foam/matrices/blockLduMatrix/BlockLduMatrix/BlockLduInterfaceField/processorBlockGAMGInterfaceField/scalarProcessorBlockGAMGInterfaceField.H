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

#ifndef scalarProcessorBlockGAMGInterfaceField_H
#define scalarProcessorBlockGAMGInterfaceField_H

#include "processorBlockGAMGInterfaceField.H"
#include "BlockGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarCoeffField.H"
#include "BlockCoeff.H"
#include "scalarBlockCoeff.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

template<>
void processorBlockGAMGInterfaceField<scalar>::initInterfaceMatrixUpdate
(
    const Field<scalar>& psiInternal,
    Field<scalar>&,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>&,
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


template<>
void processorBlockGAMGInterfaceField<scalar>::updateInterfaceMatrix
(
    const Field<scalar>& psiInternal,
    Field<scalar>& result,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    Field<scalar> pnf
    (
        procInterface_.compressedReceive<scalar>(commsType, coeffs.size())
    );

    pnf = coeffs.asScalar()*
        procInterface_.compressedReceive<scalar>(commsType, this->size())();

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
