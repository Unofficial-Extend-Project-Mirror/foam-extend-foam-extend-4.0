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

inline const enginePiston& pistonSliding::piston() const
{
    return piston_;
}


inline const thoboisSlidingValveBank& pistonSliding::valves() const
{
    return valves_;
}

inline const scalar& pistonSliding::deformSwitch() const
{
    return deformSwitch_;
}

inline const scalar& pistonSliding::delta() const
{
    return delta_;
}

inline const scalar& pistonSliding::offSet() const
{
    return offSet_;
}


inline const scalar& pistonSliding::valveTopTol() const
{
    return valveTopTol_;
}


inline const scalar& pistonSliding::pistonPosition() const
{
    return pistonPosition_;
}

inline scalar& pistonSliding::pistonPosition()
{
    return pistonPosition_;
}

inline const scalar& pistonSliding::virtualPistonPosition() const
{
    return virtualPistonPosition_;
}

inline scalar& pistonSliding::virtualPistonPosition()
{
    return virtualPistonPosition_;
}


inline const scalarField& pistonSliding::valveTopPosition() const
{
    return valveTopPosition_;
}

inline scalarField& pistonSliding::valveTopPosition()
{
    return valveTopPosition_;
}


inline const scalarField& pistonSliding::valveBottomPosition() const
{
    return valveBottomPosition_;
}

inline scalarField& pistonSliding::valveBottomPosition()
{
    return valveBottomPosition_;
}

inline const scalarField& pistonSliding::valvePistonPosition() const
{
    return valvePistonPosition_;
}

inline scalarField& pistonSliding::valvePistonPosition()
{
    return valvePistonPosition_;
}

inline const scalar& pistonSliding::deckHeight() const
{
    return deckHeight_;
}

inline scalar& pistonSliding::deckHeight()
{
    return deckHeight_;
}

inline const scalarField& pistonSliding::minValveZ() const
{
    return minValveZ_;
}

inline scalarField& pistonSliding::maxValveZ()
{
    return minValveZ_;
}

inline const scalar& pistonSliding::poppetValveTol() const
{
    return poppetValveTol_;
}

inline boolList& pistonSliding::isReallyClosed()
{
    return isReallyClosed_;
}

inline const boolList& pistonSliding::isReallyClosed() const
{
    return isReallyClosed_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
