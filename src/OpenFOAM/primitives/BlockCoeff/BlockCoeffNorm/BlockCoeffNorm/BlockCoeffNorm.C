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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    BlockCoeffNorm

\*---------------------------------------------------------------------------*/

#include "blockCoeffNorms.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeffNorm<Type>::BlockCoeffNorm
(
    const dictionary& dict
)
:
    dict_(dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<typename Foam::BlockCoeffNorm<Type> >
Foam::BlockCoeffNorm<Type>::New
(
    const dictionary& dict
)
{
    word normName(dict.lookup("norm"));

    typename dictionaryConstructorTable::iterator constructorIter =
        dictionaryConstructorTablePtr_->find(normName);

    if (constructorIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "autoPtr<BlockCoeffNorm> BlockCoeffNorm::New\n"
            "(\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown norm " << normName
            << endl << endl
            << "Valid matrix norms are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<BlockCoeffNorm<Type> >
    (
        constructorIter()
        (
            dict
        )
    );
}


// ************************************************************************* //
