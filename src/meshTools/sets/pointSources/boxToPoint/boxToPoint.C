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

#include "boxToPoint.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boxToPoint, 0);

addToRunTimeSelectionTable(topoSetSource, boxToPoint, word);

addToRunTimeSelectionTable(topoSetSource, boxToPoint, istream);

}


Foam::topoSetSource::addToUsageTable Foam::boxToPoint::usage_
(
    boxToPoint::typeName,
    "\n    Usage: boxToPoint ((minx miny minz) (maxx maxy maxz))\n\n"
    "    Select all points with coordinate within bounding box\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boxToPoint::combine(topoSet& set, const bool add) const
{
    const pointField& pts = mesh_.points();

    forAll(pts, pointI)
    {
        if (bb_.contains(pts[pointI]))
        {
            addOrDelete(set, pointI, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::boxToPoint::boxToPoint
(
    const polyMesh& mesh,
    const treeBoundBox& bb
)
:
    topoSetSource(mesh),
    bb_(bb)
{}


// Construct from dictionary
Foam::boxToPoint::boxToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    bb_(dict.lookup("box"))
{}


// Construct from Istream
Foam::boxToPoint::boxToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    bb_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boxToPoint::~boxToPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boxToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding points that are within box " << bb_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing points that are within box " << bb_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
