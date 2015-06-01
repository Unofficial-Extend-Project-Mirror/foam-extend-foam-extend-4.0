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

#ifndef tensorProcessorBlockGAMGInterfaceField_H
#define tensorProcessorBlockGAMGInterfaceField_H

#include "processorBlockGAMGInterfaceField.H"
#include "BlockGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "tensorCoeffField.H"
#include "BlockCoeff.H"
#include "tensorBlockCoeff.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

template<>
void processorBlockGAMGInterfaceField<tensor>::initInterfaceMatrixUpdate
(
    const Field<tensor>& psiInternal,
    Field<tensor>&,
    const BlockLduMatrix<tensor>&,
    const CoeffField<tensor>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    FatalErrorIn
    (
        "tensor processorBlockGAMGInterfaceField<tensor>::"
        "initInterafaceMatrix(...)"
    )   << "Not implemented" << abort(FatalError);
}

template<>
void processorBlockGAMGInterfaceField<tensor>::updateInterfaceMatrix
(
    const Field<tensor>& psiInternal,
    Field<tensor>& result,
    const BlockLduMatrix<tensor>&,
    const CoeffField<tensor>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    FatalErrorIn
    (
        "tensor processorBlockGAMGInterfaceField<tensor>::"
        "updateInterfaceMatrix(...)"
    )   << "Not implemented" << abort(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
