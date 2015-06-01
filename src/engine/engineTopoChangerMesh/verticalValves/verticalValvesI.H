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
