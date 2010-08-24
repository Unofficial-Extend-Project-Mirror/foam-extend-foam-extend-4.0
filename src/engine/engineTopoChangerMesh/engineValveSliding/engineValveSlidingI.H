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
