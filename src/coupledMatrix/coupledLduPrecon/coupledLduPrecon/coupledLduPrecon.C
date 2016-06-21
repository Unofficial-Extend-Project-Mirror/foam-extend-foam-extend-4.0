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
    coupledLduPrecon

Description
    Virtual base class for coupled matrix preconditioners

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "coupledLduPrecon.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(coupledLduPrecon, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::word Foam::coupledLduPrecon::getName
(
    const dictionary& dict
)
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


Foam::autoPtr<Foam::coupledLduPrecon> Foam::coupledLduPrecon::New
(
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
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

    dictionaryConstructorTable::iterator constructorIter =
        dictionaryConstructorTablePtr_->find(preconName);

    if (constructorIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "autoPtr<coupledLduPrecon> coupledLduPrecon::New\n"
            "(\n"
            "    const coupledLduMatrix& matrix,\n"
            "    const PtrList<FieldField<Field, scalar> >& bouCoeffs,\n"
            "    const PtrList<FieldField<Field, scalar> >& intCoeffs,\n"
            "    const lduInterfaceFieldPtrsListList& interfaces,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown matrix preconditioner " << preconName
            << endl << endl
            << "Valid matrix preconditioners are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coupledLduPrecon>
    (
        constructorIter()
        (
            matrix,
            bouCoeffs,
            intCoeffs,
            interfaces,
            controls
        )
    );
}


// ************************************************************************* //
