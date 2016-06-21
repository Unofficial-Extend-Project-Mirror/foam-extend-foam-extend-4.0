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

#include "splineEdge.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(splineEdge, 0);

    addToRunTimeSelectionTable
    (
        curvedEdge,
        splineEdge,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::splineEdge::splineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& internalPoints
)
:
    curvedEdge(points, start, end),
    CatmullRomSpline(appendEndPoints(points, start, end, internalPoints))
{}


Foam::splineEdge::splineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    CatmullRomSpline(appendEndPoints(points, start_, end_, pointField(is)))
{
    token t(is);
    is.putBack(t);

    // discard unused start/end tangents
    if (t == token::BEGIN_LIST)
    {
        vector tangent0Ignored(is);
        vector tangent1Ignored(is);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::splineEdge::~splineEdge()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::splineEdge::position(const scalar mu) const
{
    return CatmullRomSpline::position(mu);
}


Foam::scalar Foam::splineEdge::length() const
{
    return CatmullRomSpline::length();
}


// ************************************************************************* //
