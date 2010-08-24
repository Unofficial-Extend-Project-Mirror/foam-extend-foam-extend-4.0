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
