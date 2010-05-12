/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const enginePiston& accordionEngineMesh::piston() const
{
    return piston_;
}


inline const thoboisValveBank& accordionEngineMesh::valves() const
{
    return valves_;
}
                 
inline const scalar& accordionEngineMesh::deformSwitch() const
{
    return deformSwitch_;
}

inline const scalar& accordionEngineMesh::delta() const
{
    return delta_;
}

inline const scalar& accordionEngineMesh::offSet() const
{
    return offSet_;
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
