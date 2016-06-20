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
    coupledLduSmoother

Description
    Coupled LDU matrix smoother virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "coupledLduSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(coupledLduSmoother, word);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::word Foam::coupledLduSmoother::getName
(
    const dictionary& dict
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        e.dict().lookup("smoother") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::coupledLduSmoother> Foam::coupledLduSmoother::New
(
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const dictionary& dict
)
{
    word smootherName;

    // Handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        e.dict().lookup("smoother") >> smootherName;
    }
    else
    {
        e.stream() >> smootherName;
    }

    // Not (yet?) needed:
    // const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    wordConstructorTable::iterator constructorIter =
        wordConstructorTablePtr_->find(smootherName);

    if (constructorIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "(\n"
            "    const word& smootherName,\n"
            "    const coupledLduMatrix& matrix,\n"
            "    const PtrList<FieldField<Field, scalar> >& bouCoeffs,\n"
            "    const PtrList<FieldField<Field, scalar> >& intCoeffs,\n"
            "    const lduInterfaceFieldPtrsListList& interfaces\n"
            ")"
        )   << "Unknown matrix smoother " << smootherName
            << endl << endl
            << "Valid matrix smoothers are :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<coupledLduSmoother>
    (
        constructorIter()
        (
            matrix,
            bouCoeffs,
            intCoeffs,
            interfaces
        )
    );
}


// ************************************************************************* //
