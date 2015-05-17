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

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const enginePiston& accordionEngineMesh::piston() const
{
    return piston_;
}


inline const accordionValveBank& accordionEngineMesh::valves() const
{
    return valves_;
}


inline const scalar& accordionEngineMesh::pistonPosition() const
{
    return pistonPosition_;
}

inline scalar& accordionEngineMesh::pistonPosition()
{
    return pistonPosition_;
}

inline const scalar& accordionEngineMesh::virtualPistonPosition() const
{
    return virtualPistonPosition_;
}

inline scalar& accordionEngineMesh::virtualPistonPosition()
{
    return virtualPistonPosition_;
}


inline const scalar& accordionEngineMesh::deckHeight() const
{
    return deckHeight_;
}

inline scalar& accordionEngineMesh::deckHeight()
{
    return deckHeight_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
