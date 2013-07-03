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

    pnf = coeffs.asScalar() *
        procInterface_.compressedReceive<scalar>(commsType, this->size())();

    const unallocLabelList& faceCells = procInterface_.faceCells();


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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
