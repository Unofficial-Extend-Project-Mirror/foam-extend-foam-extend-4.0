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

\*---------------------------------------------------------------------------*/

#include "setToPoint.H"
#include "polyMesh.H"
#include "pointSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(setToPoint, 0);

addToRunTimeSelectionTable(topoSetSource, setToPoint, word);

addToRunTimeSelectionTable(topoSetSource, setToPoint, istream);

}


Foam::topoSetSource::addToUsageTable Foam::setToPoint::usage_
(
    setToPoint::typeName,
    "\n    Usage: setToPoint set\n\n"
    "    Select all points in the pointSet\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::setToPoint::combine(topoSet& set, const bool add) const
{
    pointSet cs
    (
        mesh_,
        setName_
    );

    const labelList pointLabels = cs.toc();

    if (pointLabels.size() > 0)
    {
        forAll (pointLabels, i)
        {
            // Only do active points
            if (pointLabels[i] < mesh_.nPoints())
            {
                addOrDelete(set, pointLabels[i], add);
            }
        }
    }
    else
    {
        WarningIn("setToPoint::combine(topoSet&, const bool)")
            << "Point set named " << setName_ << " is empty" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::setToPoint::setToPoint
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetSource(mesh),
    setName_(setName)
{}


// Construct from dictionary
Foam::setToPoint::setToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("name"))
{}


// Construct from Istream
Foam::setToPoint::setToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setToPoint::~setToPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all points of pointSet " << setName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all points of pointSet " << setName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
