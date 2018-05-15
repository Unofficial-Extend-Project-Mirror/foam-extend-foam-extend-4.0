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

#include "AMGInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AMGInterface, 0);
    defineRunTimeSelectionTable(AMGInterface, lduInterface);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::AMGInterface::combine(const AMGInterface& coarseGi)
{
    const labelList& coarseFra = coarseGi.restrictAddressing_;

    forAll (restrictAddressing_, ffi)
    {
        restrictAddressing_[ffi] = coarseFra[restrictAddressing_[ffi]];
    }

    faceCells_ = coarseGi.faceCells_;
}


Foam::tmp<Foam::labelField> Foam::AMGInterface::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return interfaceInternalField<label>(internalData);
}


Foam::tmp<Foam::scalarField> Foam::AMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    tmp<scalarField> tcoarseCoeffs(new scalarField(size(), 0.0));
    scalarField& coarseCoeffs = tcoarseCoeffs();

    // Added weights to account for non-integral matching
    forAll (restrictAddressing_, ffi)
    {
        coarseCoeffs[restrictAddressing_[ffi]] +=
            restrictWeights_[ffi]*fineCoeffs[fineAddressing_[ffi]];
    }

    return tcoarseCoeffs;
}


// ************************************************************************* //
