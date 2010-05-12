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

Description

\*---------------------------------------------------------------------------*/

#include "cellToFace.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(cellToFace, 0);

addToRunTimeSelectionTable(topoSetSource, cellToFace, word);

addToRunTimeSelectionTable(topoSetSource, cellToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::cellToFace::usage_
(
    cellToFace::typeName,
    "\n    Usage: cellToFace <cellSet> all|both\n\n"
    "    Select -all : all faces of cells in the cellSet\n"
    "           -both: faces where both neighbours are in the cellSet\n\n"
);

template<>
const char* Foam::NamedEnum<Foam::cellToFace::cellAction, 2>::names[] =
{
    "all",
    "both"
};

const Foam::NamedEnum<Foam::cellToFace::cellAction, 2>
    Foam::cellToFace::cellActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToFace::combine(topoSet& set, const bool add) const
{
    // Load the set
    if (!exists(mesh_.time().path()/topoSet::localPath(mesh_, setName_)))
    {
        SeriousError<< "Cannot load set "
            << setName_ << endl;
    }
    
    cellSet loadedSet(mesh_, setName_);

    if (option_ == ALL)
    {
        // Add all faces from cell
        for
        (
            cellSet::const_iterator iter = loadedSet.begin();
            iter != loadedSet.end();
            ++iter
        )
        {
            label cellI = iter.key();

            const labelList& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                addOrDelete(set, cFaces[cFaceI], add);
            }
        }
    }
    else if (option_ == BOTH)
    {
        // Add all faces whose both neighbours are in set.

        // Count number of cells using face.
        Map<label> numCells(loadedSet.size());

        for
        (
            cellSet::const_iterator iter = loadedSet.begin();
            iter != loadedSet.end();
            ++iter
        )
        {
            label cellI = iter.key();

            const labelList& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];

                Map<label>::iterator fndFace = numCells.find(faceI);

                if (fndFace == numCells.end())
                {
                    numCells.insert(faceI, 1);
                }
                else
                {
                    fndFace()++;
                }
            }
        }

        // Include faces that are referenced twice
        for
        (
            Map<label>::const_iterator iter = numCells.begin();
            iter != numCells.end();
            ++iter
        )
        {
            if (iter() == 2)
            {
                addOrDelete(set, iter.key(), add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from componenta
Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const word& setName,
    const cellAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


// Construct from dictionary
Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const dictionary& dict          
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(cellActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(cellActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellToFace::~cellToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Pout<< "    Adding faces according to cellSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Pout<< "    Removing faces according to cellSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
