/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coordinateModification.H"
#include "dictionary.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<Foam::coordinateModification> Foam::coordinateModification::New
(
    const word& name,
    const dictionary& dict
)
{
    if( debug )
    {
        Info<< "coordinateModification::New(const word&, const dictionary&) : "
            << "constructing coordinateModification"
            << endl;
    }

    // default type is self
    word cmType(typeName_());
    if( dict.found("type") )
    {
        dict.lookup("type") >> cmType;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(cmType);

    if( cstrIter == dictionaryConstructorTablePtr_->end() )
    {
        FatalIOErrorIn
        (
            "coordinateModification::New(const word&, const dictionary&)",
            dict
        )   << "Unknown coordinateModification type " << cmType << nl << nl
            << "Valid coordinateModification types are :" << nl
            << "[default: " << typeName_() << "]"
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateModification>(cstrIter()(name, dict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
