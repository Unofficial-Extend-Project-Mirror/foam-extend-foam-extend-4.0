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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::curvedEdge::start() const
{
    return start_;
}


inline Foam::label Foam::curvedEdge::end() const
{
    return end_;
}


inline int Foam::curvedEdge::compare(const label start, const label end) const
{
    if (start_ == start && end_ == end)
    {
        return 1;
    }
    else if (start_ == end && end_ == start)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}


inline int Foam::curvedEdge::compare(const curvedEdge& e) const
{
    return Foam::curvedEdge::compare(e.start(), e.end());
}


inline int Foam::curvedEdge::compare(const edge& e) const
{
    return Foam::curvedEdge::compare(e.start(), e.end());
}


// ************************************************************************* //
