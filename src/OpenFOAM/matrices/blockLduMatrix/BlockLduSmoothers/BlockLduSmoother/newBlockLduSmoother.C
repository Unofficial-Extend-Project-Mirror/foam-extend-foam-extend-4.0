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
    BlockLduSmoother

Description
    Block LDU matrix smoother virtual base class

\*----------------------------------------------------------------------------*/

#include "BlockLduSmoother.H"

template<class Type>
class BlockNoSmoother;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduSmoother<Type> > Foam::BlockLduSmoother<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const FieldField<CoeffField, Type>& boundaryCoeffs,
    const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces,
    const dictionary& dict
)
{
    word smootherName(dict.lookup("smoother"));

    typename dictionaryConstructorTable::iterator constructorIter =
        dictionaryConstructorTablePtr_->find(smootherName);

    if (constructorIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "autoPtr<BlockLduSmoother> BlockLduSmoother::New\n"
            "(\n"
            "    const BlockLduMatrix<Type>& matrix,\n"
            "    const FieldField<CoeffField, Type>& boundaryCoeffs,\n"
            "    const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown matrix smoother " << smootherName
            << endl << endl
            << "Valid matrix smoothers are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<BlockLduSmoother<Type> >
    (
        constructorIter()
        (
            matrix,
            boundaryCoeffs,
            interfaces,
            dict
        )
    );
}


// ************************************************************************* //
