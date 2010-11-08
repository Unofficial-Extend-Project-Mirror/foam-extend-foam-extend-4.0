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

\*---------------------------------------------------------------------------*/

#include "BlockLduPrecon.H"
#include "blockNoPrecons.H"

template<class Type>
class BlockNoPrecon;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduPrecon<Type> > Foam::BlockLduPrecon<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
{
    word preconName;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> preconName;
    }
    else
    {
        e.stream() >> preconName;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (matrix.diagonal())
    {
        // No preconditioning for the diagonal matrix
        return autoPtr<BlockLduPrecon<Type> >
        (
            new BlockNoPrecon<Type>
            (
                matrix,
                controls
            )
        );
    }
    else
    {
        typename dictionaryConstructorTable::iterator constructorIter =
            dictionaryConstructorTablePtr_->find(preconName);

        if (constructorIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "autoPtr<BlockLduPrecon> BlockLduPrecon::New\n"
                "(\n"
                "    const BlockLduMatrix<Type>& matrix,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown matrix preconditioner " << preconName
                << endl << endl
                << "Valid matrix preconditioners are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<BlockLduPrecon<Type> >
        (
            constructorIter()
            (
                matrix,
                controls
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::BlockLduPrecon<Type>::getName(const dictionary& dict)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


// ************************************************************************* //
