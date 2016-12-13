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

#include "pointToCell.H"
#include "polyMesh.H"
#include "pointSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pointToCell, 0);

addToRunTimeSelectionTable(topoSetSource, pointToCell, word);

addToRunTimeSelectionTable(topoSetSource, pointToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::pointToCell::usage_
(
    pointToCell::typeName,
    "\n    Usage: pointToCell <pointSet> any|all\n\n"
    "    Select cells with\n"
    "    -any point in the pointSet\n"
    "    -all points in the pointSet\n\n"
);

template<>
const char* Foam::NamedEnum<Foam::pointToCell::pointAction, 2>::names[] =
{
    "any",
    "all"
};


const Foam::NamedEnum<Foam::pointToCell::pointAction, 2>
    Foam::pointToCell::pointActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointToCell::combine(topoSet& set, const bool add) const
{
    // Load the set
    pointSet loadedSet(mesh_, setName_);


    if (option_ == ANY)
    {
        // Add cells with any point in loadedSet
        for
        (
            pointSet::const_iterator iter = loadedSet.begin();
            iter != loadedSet.end();
            ++iter
        )
        {
            label pointI = iter.key();

            const labelList& pCells = mesh_.pointCells()[pointI];

            forAll(pCells, pCellI)
            {
                addOrDelete(set, pCells[pCellI], add);
            }
        }
    }
    else if (option_ == ALL)
    {
        // Add all cells whose points are all in set.

        Map<label> numPoints(loadedSet.size());

        forAllConstIter(pointSet, loadedSet, iter)
        {
            label pointI = iter.key();

            const labelList& pCells = mesh_.pointCells()[pointI];

            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];

                Map<label>::iterator fndCell = numPoints.find(cellI);

                if (fndCell == numPoints.end())
                {
                    numPoints.insert(cellI, 1);
                }
                else
                {
                    fndCell()++;
                }
            }
        }


        // Include cells that are referenced as many times as there are points
        // in cell -> all points of cell
        for
        (
            Map<label>::const_iterator iter = numPoints.begin();
            iter != numPoints.end();
            ++iter
        )
        {
            label cellI = iter.key();

            if (iter() == mesh_.cellPoints()[cellI].size())
            {
                addOrDelete(set, cellI, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


// Construct from dictionary
Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(pointActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointToCell::~pointToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
