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

Class
    lduPrecon

Description
    Virtual base class for LDU preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "amgPolicy.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable
    (
        amgPolicy,
        matrix
    );

} // End namespace Foam


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::amgPolicy> Foam::amgPolicy::New
(
    const word& policyType,
    const lduMatrix& matrix,
    const label groupSize,
    const label nCoarseCells
)
{
    matrixConstructorTable::iterator constructorIter =
        matrixConstructorTablePtr_->find(policyType);

    if (constructorIter == matrixConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "autoPtr<amgPolicy> amgPolicy::New\n"
            "(\n"
            "    const word& policyType,\n"
            "    const lduMatrix& matrix,\n"
            "    const label groupSize\n"
            "    const label nCoarseCells\n"
            ")"
        )   << "Unknown AMG policy " << policyType
            << endl << endl
            << "Valid AMG policies are :" << endl
            << matrixConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<amgPolicy>
    (
        constructorIter()
        (
            matrix,
            groupSize,
            nCoarseCells
        )
    );
}


// ************************************************************************* //
