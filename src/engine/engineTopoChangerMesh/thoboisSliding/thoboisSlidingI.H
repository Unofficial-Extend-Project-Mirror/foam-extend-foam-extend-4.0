/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2010 Tommaso Lucchini
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
