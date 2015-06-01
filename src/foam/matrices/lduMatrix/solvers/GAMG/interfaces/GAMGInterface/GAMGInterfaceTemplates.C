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

#include "GAMGInterface.H"
#include "scalarCoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::GAMGInterface::interfaceInternalField
(
    const UList<Type>& iF
) const
{
    tmp<Field<Type> > tresult(new Field<Type>(size()));
    Field<Type>& result = tresult();

    forAll (result, elemI)
    {
        result[elemI] = iF[faceCells_[elemI]];
    }

    return tresult;
}

template<class Type>
Foam::tmp<Foam::CoeffField<Type> > Foam::GAMGInterface::agglomerateBlockCoeffs
(
    const Foam::CoeffField<Type>& fineCoeffs
) const
{
    tmp<CoeffField<Type> > tcoarseCoeffs(new CoeffField<Type>(size()));
    CoeffField<Type>& coarseCoeffs = tcoarseCoeffs();

    typedef CoeffField<Type> TypeCoeffField;

    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    if (fineCoeffs.activeType() == blockCoeffBase::SQUARE)
    {
        squareTypeField& activeCoarseCoeffs = coarseCoeffs.asSquare();
        const squareTypeField& activeFineCoeffs = fineCoeffs.asSquare();

        activeCoarseCoeffs *= 0.0;

        // Added weights to account for non-integral matching
        forAll(restrictAddressing_, ffi)
        {
            activeCoarseCoeffs[restrictAddressing_[ffi]] +=
                restrictWeights_[ffi]*activeFineCoeffs[fineAddressing_[ffi]];
        }
    }
    else if (fineCoeffs.activeType() == blockCoeffBase::LINEAR)
    {
        linearTypeField& activeCoarseCoeffs = coarseCoeffs.asLinear();
        const linearTypeField& activeFineCoeffs = fineCoeffs.asLinear();

        activeCoarseCoeffs *= 0.0;

        // Added weights to account for non-integral matching
        forAll(restrictAddressing_, ffi)
        {
            activeCoarseCoeffs[restrictAddressing_[ffi]] +=
                restrictWeights_[ffi]*activeFineCoeffs[fineAddressing_[ffi]];
        }
    }

    return tcoarseCoeffs;
}

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
