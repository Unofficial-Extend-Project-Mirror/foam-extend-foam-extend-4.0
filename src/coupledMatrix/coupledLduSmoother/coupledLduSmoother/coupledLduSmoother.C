/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    coupledLduSmoother

Description
    Coupled LDU matrix smoother virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "coupledLduSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(coupledLduSmoother, word);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coupledLduSmoother> Foam::coupledLduSmoother::New
(
    const word& smootherName,
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces
)
{
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
            << wordConstructorTablePtr_->toc()
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
