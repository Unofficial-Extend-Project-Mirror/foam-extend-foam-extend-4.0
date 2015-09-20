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

inline const enginePiston& thoboisSliding::piston() const
{
    return piston_;
}


inline const thoboisSlidingValveBank& thoboisSliding::valves() const
{
    return valves_;
}

inline const scalar& thoboisSliding::deformSwitch() const
{
    return deformSwitch_;
}

inline const scalar& thoboisSliding::delta() const
{
    return delta_;
}

inline const scalar& thoboisSliding::offSet() const
{
    return offSet_;
}


inline const scalar& thoboisSliding::valveTopTol() const
{
    return valveTopTol_;
}


inline const scalar& thoboisSliding::pistonPosition() const
{
    return pistonPosition_;
}

inline scalar& thoboisSliding::pistonPosition()
{
    return pistonPosition_;
}

inline const scalar& thoboisSliding::virtualPistonPosition() const
{
    return virtualPistonPosition_;
}

inline scalar& thoboisSliding::virtualPistonPosition()
{
    return virtualPistonPosition_;
}


inline const scalarField& thoboisSliding::valveTopPosition() const
{
    return valveTopPosition_;
}

inline scalarField& thoboisSliding::valveTopPosition()
{
    return valveTopPosition_;
}


inline const scalarField& thoboisSliding::valveBottomPosition() const
{
    return valveBottomPosition_;
}

inline scalarField& thoboisSliding::valveBottomPosition()
{
    return valveBottomPosition_;
}

inline const scalarField& thoboisSliding::valvePistonPosition() const
{
    return valvePistonPosition_;
}

inline scalarField& thoboisSliding::valvePistonPosition()
{
    return valvePistonPosition_;
}

inline const scalar& thoboisSliding::deckHeight() const
{
    return deckHeight_;
}

inline scalar& thoboisSliding::deckHeight()
{
    return deckHeight_;
}

inline const scalarField& thoboisSliding::minValveZ() const
{
    return minValveZ_;
}

inline scalarField& thoboisSliding::maxValveZ()
{
    return minValveZ_;
}

inline const scalar& thoboisSliding::poppetValveTol() const
{
    return poppetValveTol_;
}

inline boolList& thoboisSliding::isReallyClosed()
{
    return isReallyClosed_;
}

inline const boolList& thoboisSliding::isReallyClosed() const
{
    return isReallyClosed_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
