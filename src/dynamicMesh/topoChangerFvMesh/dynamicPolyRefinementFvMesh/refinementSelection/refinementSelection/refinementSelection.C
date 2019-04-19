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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "refinementSelection.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementSelection, 0);
defineRunTimeSelectionTable(refinementSelection, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementSelection::refinementSelection
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    coeffDict_(dict.subDict("refinementSelection"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::refinementSelection> Foam::refinementSelection::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Get subdictionary from the dictionary
    const dictionary coeffDict(dict.subDict("refinementSelection"));

    // Get the name of the desired refinement selection algorithm
    const word refinementSelectionTypeName(coeffDict.lookup("type"));
    Info<< "Creating refinementSelection " << refinementSelectionTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(refinementSelectionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
            "refinementSelection::refinementSelection::New\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
            )   << "Unknown refinementSelection type "
                << refinementSelectionTypeName << endl << endl
                << "Valid refinementSelection types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalError);
    }

    return autoPtr<refinementSelection>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::refinementSelection::~refinementSelection()
{}


// ************************************************************************* //
