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
    verticalValves

\*---------------------------------------------------------------------------*/

namespace Foam
{

inline const simpleEnginePiston& verticalValves::piston() const
{
    return piston_;
}


inline const verticalValveBank& verticalValves::valves() const
{
    return valves_;
}
                 
inline const scalar& verticalValves::deformSwitch() const
{
    return deformSwitch_;
}

inline const scalar& verticalValves::delta() const
{
    return delta_;
}

inline const scalar& verticalValves::offSet() const
{
    return offSet_;
}


inline const scalar& verticalValves::valveTopTol() const
{
    return valveTopTol_;
}


inline const scalar& verticalValves::pistonPosition() const
{
    return pistonPosition_;
}

inline scalar& verticalValves::pistonPosition()
{
    return pistonPosition_;
}

inline const scalar& verticalValves::virtualPistonPosition() const
{
    return virtualPistonPosition_;
}

inline scalar& verticalValves::virtualPistonPosition()
{
    return virtualPistonPosition_;
}


inline const scalarField& verticalValves::valveTopPosition() const
{
    return valveTopPosition_;
}

inline scalarField& verticalValves::valveTopPosition()
{
    return valveTopPosition_;
}


inline const scalarField& verticalValves::valveBottomPosition() const
{
    return valveBottomPosition_;
}

inline scalarField& verticalValves::valveBottomPosition()
{
    return valveBottomPosition_;
}

inline const scalarField& verticalValves::valvePistonPosition() const
{
    return valvePistonPosition_;
}

inline scalarField& verticalValves::valvePistonPosition()
{
    return valvePistonPosition_;
}

inline const scalar& verticalValves::deckHeight() const
{
    return deckHeight_;
}

inline scalar& verticalValves::deckHeight()
{
    return deckHeight_;
}

inline const scalarField& verticalValves::minValveZ() const
{
    return minValveZ_;
}

inline scalarField& verticalValves::maxValveZ()
{
    return minValveZ_;
}

inline const scalar& verticalValves::poppetValveTol() const
{
    return poppetValveTol_;
}

inline boolList& verticalValves::isReallyClosed()
{
    return isReallyClosed_;
}

inline const boolList& verticalValves::isReallyClosed() const
{
    return isReallyClosed_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
