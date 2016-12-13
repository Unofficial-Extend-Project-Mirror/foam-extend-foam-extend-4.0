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

#include "BlockAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockAMGInterfaceField<Type> >
Foam::BlockAMGInterfaceField<Type>::New
(
    const AMGInterface& AMGCp,
    const BlockLduInterfaceField<Type>& fineInterface
)
{
    word coupleType(fineInterface.interfaceFieldType());

    typename lduInterfaceConstructorTable::iterator cstrIter =
        lduInterfaceConstructorTablePtr_->find(coupleType);

    if (cstrIter == lduInterfaceConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "BlockAMGInterfaceField::New"
            "(const AMGInterface& AMGCp, "
            "const BlockLduInterfaceField<Type>& fineInterface)"
        )   << "Unknown BlockAMGInterfaceField type " << coupleType << ".\n"
            << "Valid BlockAMGInterfaceField types are :"
            << lduInterfaceConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<BlockAMGInterfaceField<Type> >
    (
        cstrIter()
        (
            AMGCp,
            fineInterface
        )
    );
}


// ************************************************************************* //
