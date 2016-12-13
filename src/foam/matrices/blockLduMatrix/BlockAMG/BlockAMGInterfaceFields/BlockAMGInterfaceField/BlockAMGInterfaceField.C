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

\*---------------------------------------------------------------------------*/

#include "BlockAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::CoeffField<Type> >
Foam::BlockAMGInterfaceField<Type>::agglomerateBlockCoeffs
(
    const Foam::CoeffField<Type>& fineCoeffs
) const
{
    tmp<CoeffField<Type> > tcoarseCoeffs
    (
        new CoeffField<Type>(interface_.size())
    );
    CoeffField<Type>& coarseCoeffs = tcoarseCoeffs();

    typedef CoeffField<Type> TypeCoeffField;

    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    // Get addressing from the fine interface
    const labelField& fineAddressing = interface_.fineAddressing();
    const labelField& restrictAddressing = interface_.restrictAddressing();
    const scalarField& restrictWeights = interface_.restrictWeights();

    // Added weights to account for non-integral matching
    if (fineCoeffs.activeType() == blockCoeffBase::SQUARE)
    {
        squareTypeField& activeCoarseCoeffs = coarseCoeffs.asSquare();
        const squareTypeField& activeFineCoeffs = fineCoeffs.asSquare();

        activeCoarseCoeffs *= 0.0;

        // Added weights to account for non-integral matching
        forAll (restrictAddressing, ffi)
        {
            activeCoarseCoeffs[restrictAddressing[ffi]] +=
                restrictWeights[ffi]*activeFineCoeffs[fineAddressing[ffi]];
        }
    }
    else if (fineCoeffs.activeType() == blockCoeffBase::LINEAR)
    {
        linearTypeField& activeCoarseCoeffs = coarseCoeffs.asLinear();
        const linearTypeField& activeFineCoeffs = fineCoeffs.asLinear();

        activeCoarseCoeffs *= 0.0;

        // Added weights to account for non-integral matching
        forAll (restrictAddressing, ffi)
        {
            activeCoarseCoeffs[restrictAddressing[ffi]] +=
                restrictWeights[ffi]*activeFineCoeffs[fineAddressing[ffi]];
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::CoeffField<Type> >\n"
            "Foam::BlockAMGInterfaceField<Type>::agglomerateBlockCoeffs\n"
            "(\n"
            "    const Foam::CoeffField<Type>& fineCoeffs\n"
            ") const"
        )   << "Scalar type agglomeration currently not handled"
            << abort(FatalError);
    }

    return tcoarseCoeffs;
}


// ************************************************************************* //
