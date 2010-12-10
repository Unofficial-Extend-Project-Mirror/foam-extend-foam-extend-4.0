/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
            << matrixConstructorTablePtr_->toc()
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
