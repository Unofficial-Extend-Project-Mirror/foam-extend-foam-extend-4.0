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

#include "setToFace.H"
#include "polyMesh.H"
#include "faceSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(setToFace, 0);

addToRunTimeSelectionTable(topoSetSource, setToFace, word);

addToRunTimeSelectionTable(topoSetSource, setToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::setToFace::usage_
(
    setToFace::typeName,
    "\n    Usage: setToFace set\n\n"
    "    Select all faces in the faceSet\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::setToFace::combine(topoSet& set, const bool add) const
{
    faceSet cs
    (
        mesh_,
        setName_
    );

    const labelList faceLabels = cs.toc();

    if (faceLabels.size() > 0)
    {
        forAll (faceLabels, i)
        {
            // Only do active faces
            if (faceLabels[i] < mesh_.nFaces())
            {
                addOrDelete(set, faceLabels[i], add);
            }
        }
    }
    else
    {
        WarningIn("setToFace::combine(topoSet&, const bool)")
            << "Face set named " << setName_ << " is empty" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::setToFace::setToFace
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetSource(mesh),
    setName_(setName)
{}


// Construct from dictionary
Foam::setToFace::setToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("name"))
{}


// Construct from Istream
Foam::setToFace::setToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setToFace::~setToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all faces of faceSet " << setName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all faces of faceSet " << setName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
