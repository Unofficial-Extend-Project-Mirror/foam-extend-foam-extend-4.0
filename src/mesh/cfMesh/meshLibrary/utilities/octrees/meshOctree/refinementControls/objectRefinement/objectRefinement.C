/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "IOstream.H"
#include "objectRefinement.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(objectRefinement, 0);
defineRunTimeSelectionTable(objectRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectRefinement::objectRefinement()
:
    name_(),
    cellSize_(),
    refThickness_(0.0)
{}


objectRefinement::objectRefinement
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    cellSize_(readScalar(dict.lookup("cellSize"))),
    refThickness_(0.0)
{
    if( dict.found("refinementThickness") )
        refThickness_ = readScalar(dict.lookup("refinementThickness"));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

objectRefinement::~objectRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const objectRefinement& obr)
{
    os << obr.name() << nl;
    obr.writeDict(os, true);
    return os;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
