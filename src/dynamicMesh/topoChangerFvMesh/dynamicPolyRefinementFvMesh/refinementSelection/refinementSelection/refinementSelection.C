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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "refinementSelection.H"
#include "dynamicPolyRefinementFvMesh.H"
#include "polyhedralRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementSelection, 0);
defineRunTimeSelectionTable(refinementSelection, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementSelection::refinementSelection
(
    const dynamicPolyRefinementFvMesh& dynamicRefMesh,
    const polyhedralRefinement& pRef,
    const dictionary& dict
)
:
    dynamicRefMesh_(dynamicRefMesh),
    pRef_(pRef),
    coeffDict_
    (
        dict.subDict("refinementSelection")
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::refinementSelection> Foam::refinementSelection::New
(
    const dynamicPolyRefinementFvMesh& dynamicRefMesh,
    const polyhedralRefinement& pRef,
    const dictionary& dict
)
{
    // Get subdictionary from the dictionary
    const dictionary coeffDict(dict.subDict("refinementSelection"));

    // Get the name of the desired refinement selection algorithm
    const word refinementSelectionTypeName(coeffDict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(refinementSelectionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
            "refinementSelection::refinementSelection::New\n"
            "(\n"
            "    const dynamicPolyRefinementFvMesh& dynamicRefMesh,\n"
            "    const polyhedralRefinement& pRef,\n"
            "    const dictionary& dict\n"
            ")"
            )   << "Unknown refinementSelection type "
                << refinementSelectionTypeName << endl << endl
                << "Valid refinementSelection types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalError);
    }

    return
        autoPtr<refinementSelection>(cstrIter()(dynamicRefMesh, pRef, dict));
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::refinementSelection::~refinementSelection()
{}


// ************************************************************************* //
