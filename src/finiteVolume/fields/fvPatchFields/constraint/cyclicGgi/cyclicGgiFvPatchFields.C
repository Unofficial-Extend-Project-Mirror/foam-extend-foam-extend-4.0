/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Generalized grid interface (cyclicGgi) patch field, providing coupling
    between arbitrary patches

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor:
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "scalarCoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
void cyclicGgiFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    const Field<scalar>& psiInternal,
    Field<scalar>& result,
    const BlockLduMatrix<scalar>& m,
    const CoeffField<scalar>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = cyclicGgiPatch_.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    // Note: scalar interpolate does not get a transform, so this is safe
    // HJ, 12/Jan/2009
    scalarField pnf = cyclicGgiPatch_.interpolate(sField);

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = cyclicGgiPatch_.faceCells();

    forAll(fc, elemI)
    {
        result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
    };
}


template<>
void cyclicGgiFvPatchField<vector>::initInterfaceMatrixUpdate
(
    const Field<vector>& psiInternal,
    Field<vector>& result,
    const BlockLduMatrix<vector>& m,
    const CoeffField<vector>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = cyclicGgiPatch_.shadow().faceCells();

    Field<vector> sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    // Transformation is handled in interpolation.  HJ, 7/Jan/2009
    Field<vector> pnf = cyclicGgiPatch_.interpolate(sField);

    if (coeffs.activeType() == blockCoeffBase::SCALAR)
    {
        pnf = coeffs.asScalar() * pnf;
    }
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)
    {
        pnf = cmptMultiply(coeffs.asLinear(), pnf);
    }
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)
    {
        pnf = coeffs.asSquare() & pnf;
    }

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = cyclicGgiPatch_.faceCells();

    // Multiply the field by coefficients and add into the result
    forAll(fc, elemI)
    {
        result[fc[elemI]] -= pnf[elemI];
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(cyclicGgi);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
