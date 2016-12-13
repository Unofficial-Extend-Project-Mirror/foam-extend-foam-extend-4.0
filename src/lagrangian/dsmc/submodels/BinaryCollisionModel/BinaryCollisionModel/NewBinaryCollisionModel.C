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

#include "BinaryCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::BinaryCollisionModel<CloudType> >
Foam::BinaryCollisionModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word BinaryCollisionModelType(dict.lookup("BinaryCollisionModel"));

    Info<< "Selecting BinaryCollisionModel "
        << BinaryCollisionModelType
        << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(BinaryCollisionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "BinaryCollisionModel<CloudType>::New"
            "(const dictionary&, CloudType&)"
        )
            << "Unknown BinaryCollisionModelType type "
            << BinaryCollisionModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid BinaryCollisionModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<BinaryCollisionModel<CloudType> >
    (
        cstrIter()(dict, owner)
    );
}


// ************************************************************************* //
