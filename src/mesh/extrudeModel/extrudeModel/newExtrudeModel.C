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

#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extrudeModel> Foam::extrudeModel::New
(
    const dictionary& dict
)
{
    word extrudeModelType(dict.lookup("extrudeModel"));

    Info<< "Selecting extrudeModel " << extrudeModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(extrudeModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("extrudeModel::New(const dictionary&)")
            << "Unknown extrudeModelType type "
            << extrudeModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid extrudeModel types are :" << nl
            << dictionaryConstructorTablePtr_->toc() << nl
            << exit(FatalError);
    }

    return autoPtr<extrudeModel>(cstrIter()(dict));
}


// ************************************************************************* //

