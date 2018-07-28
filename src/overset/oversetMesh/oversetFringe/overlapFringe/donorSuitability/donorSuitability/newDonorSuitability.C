/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "error.H"
#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::donorSuitability::donorSuitability>
Foam::donorSuitability::donorSuitability::New
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
{
    dictionary coeffDict(dict.subDict("donorSuitability"));

    word donorSuitabilityTypeName(coeffDict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(donorSuitabilityTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "donorSuitability::donorSuitability::New\n"
            "(\n"
            "    const oversetFringe& oversetFringeAlgorithm,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unknown donorSuitability type " << donorSuitabilityTypeName
            << endl << endl
            << "Valid donorSuitability types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<donorSuitability>(cstrIter()(oversetFringeAlgorithm, dict));
}


// ************************************************************************* //
