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

#include "error.H"
#include "polyLineEdge.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyLineEdge, 0);
    addToRunTimeSelectionTable(curvedEdge, polyLineEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyLineEdge::polyLineEdge
(
    const pointField& ps,
    const label start,
    const label end,
    const pointField& otherPoints
)
:
    curvedEdge(ps, start, end),
    polyLine(appendEndPoints(ps, start_, end_, otherPoints))
{}


Foam::polyLineEdge::polyLineEdge(const pointField& ps, Istream& is)
:
    curvedEdge(ps, is),
    polyLine(appendEndPoints(ps, start_, end_, pointField(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyLineEdge::~polyLineEdge()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::polyLineEdge::position(const scalar lambda) const
{
    return polyLine::position(lambda);
}


Foam::scalar Foam::polyLineEdge::length() const
{
    return polyLine::lineLength_;
}


// ************************************************************************* //
