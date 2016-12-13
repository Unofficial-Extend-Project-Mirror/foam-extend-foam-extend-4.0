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

#include "regionCoupleAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupleAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        AMGInterfaceField,
        regionCoupleAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupleAMGInterfaceField::regionCoupleAMGInterfaceField
(
    const AMGInterface& AMGCp,
    const lduInterfaceField& fineInterface
)
:
    ggiAMGInterfaceField(AMGCp, fineInterface),
    regionCoupleInterface_(refCast<const regionCoupleAMGInterface>(AMGCp))
{}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::regionCoupleAMGInterfaceField::~regionCoupleAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupleAMGInterfaceField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // This must have a reduce in it.  HJ, 15/May/2009
    if (regionCoupleInterface_.coupled())
    {
        ggiAMGInterfaceField::initInterfaceMatrixUpdate
        (
            psiInternal,
            result,
            m,
            coeffs,
            cmpt,
            commsType,
            switchToLhs
        );
    }
}


void Foam::regionCoupleAMGInterfaceField::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Get expanded data to zone size.  No global reduce allowed
    // HJ, 15/May/2009
    if (regionCoupleInterface_.coupled())
    {
        ggiAMGInterfaceField::updateInterfaceMatrix
        (
            psiInternal,
            result,
            m,
            coeffs,
            cmpt,
            commsType,
            switchToLhs
        );
    }
}


// ************************************************************************* //
