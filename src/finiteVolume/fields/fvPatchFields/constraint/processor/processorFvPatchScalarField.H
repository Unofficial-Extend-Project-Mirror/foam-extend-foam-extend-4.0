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

#ifndef processorFvPatchScalarField_H
#define processorFvPatchScalarField_H

#include "processorFvPatchField.H"
#include "scalarCoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void processorFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const;


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const;


template<>
void processorFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const;


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
