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
    BlockAmgPolicy

\*---------------------------------------------------------------------------*/

#include "BlockAmgPolicy.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockAmgPolicy<Type> > Foam::BlockAmgPolicy<Type>::New
(
    const word& policyType,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const label groupSize,
    const label nCoarseCells
)
{
    typename matrixConstructorTable::iterator constructorIter =
        matrixConstructorTablePtr_->find(policyType);

    if (constructorIter == matrixConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "autoPtr<BlockAmgPolicy<Type> > BlockAmgPolicy<Type>::New\n"
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

    return autoPtr<BlockAmgPolicy<Type> >
    (
        constructorIter()
        (
            matrix,
            dict,
            groupSize,
            nCoarseCells
        )
    );
}


// ************************************************************************* //
