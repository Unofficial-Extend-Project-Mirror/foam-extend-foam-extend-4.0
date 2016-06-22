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

#include "setToCell.H"
#include "polyMesh.H"
#include "cellSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(setToCell, 0);

addToRunTimeSelectionTable(topoSetSource, setToCell, word);

addToRunTimeSelectionTable(topoSetSource, setToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::setToCell::usage_
(
    setToCell::typeName,
    "\n    Usage: setToCell set\n\n"
    "    Select all cells in the cellSet\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::setToCell::combine(topoSet& set, const bool add) const
{
    cellSet cs
    (
        mesh_,
        setName_
    );

    const labelList cellLabels = cs.toc();

    if (cellLabels.size() > 0)
    {
        forAll (cellLabels, i)
        {
            // Only do active cells
            if (cellLabels[i] < mesh_.nCells())
            {
                addOrDelete(set, cellLabels[i], add);
            }
        }
    }
    else
    {
        WarningIn("setToCell::combine(topoSet&, const bool)")
            << "Cell set named " << setName_ << " is empty" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::setToCell::setToCell
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetSource(mesh),
    setName_(setName)
{}


// Construct from dictionary
Foam::setToCell::setToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("name"))
{}


// Construct from Istream
Foam::setToCell::setToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setToCell::~setToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells of cellSet " << setName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of cellSet " << setName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
