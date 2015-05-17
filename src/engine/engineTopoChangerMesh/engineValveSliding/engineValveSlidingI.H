/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Class
    verticalValvesGambit

\*---------------------------------------------------------------------------*/

namespace Foam
{

inline const enginePiston& engineValveSliding::piston() const
{
    return piston_;
}


inline const thoboisSlidingValveBank& engineValveSliding::valves() const
{
    return valves_;
}

inline const scalar& engineValveSliding::deformSwitch() const
{
    return deformSwitch_;
}


inline const scalar& engineValveSliding::valveTopTol() const
{
    return valveTopTol_;
}


inline const scalar& engineValveSliding::pistonPosition() const
{
    return pistonPosition_;
}

inline scalar& engineValveSliding::pistonPosition()
{
    return pistonPosition_;
}

inline const scalar& engineValveSliding::deckHeight() const
{
    return deckHeight_;
}

inline scalar& engineValveSliding::deckHeight()
{
    return deckHeight_;
}

inline const scalarField& engineValveSliding::minValveZ() const
{
    return minValveZ_;
}

inline scalarField& engineValveSliding::maxValveZ()
{
    return minValveZ_;
}

inline const scalar& engineValveSliding::poppetValveTol() const
{
    return poppetValveTol_;
}

inline boolList& engineValveSliding::isReallyClosed()
{
    return isReallyClosed_;
}

inline const boolList& engineValveSliding::isReallyClosed() const
{
    return isReallyClosed_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
