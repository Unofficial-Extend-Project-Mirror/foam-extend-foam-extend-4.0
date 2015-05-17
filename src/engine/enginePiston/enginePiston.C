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

\*---------------------------------------------------------------------------*/

#include "enginePiston.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::enginePiston::enginePiston
(
    const polyMesh& mesh,
    const word& pistonPatchName,
    const autoPtr<coordinateSystem>& pistonCS,
    const scalar minLayer,
    const scalar maxLayer,
    const word& pistonPointSetName,
    const word& pistonFaceSetName,
    const word& pistonCellSetName,
    const word& bowlInPistonPatchName,
    const word& bowlInCylinderPatchName
)
:
    simpleEnginePiston
    (
        mesh,
        pistonPatchName,
        pistonCS,
        minLayer,
        maxLayer
    ),
    pistonPointSetName_(pistonPointSetName),
    pistonFaceSetName_(pistonFaceSetName),
    pistonCellSetName_(pistonPointSetName),
    bowlInPistonPatchID_(bowlInPistonPatchName, mesh.boundaryMesh()),
    bowlInCylinderPatchID_(bowlInCylinderPatchName, mesh.boundaryMesh())
{}


// Construct from dictionary
Foam::enginePiston::enginePiston
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    simpleEnginePiston(mesh, dict),
    pistonPointSetName_(dict.lookup("pistonPointSetName")),
    pistonFaceSetName_(dict.lookup("pistonFaceSetName")),
    pistonCellSetName_(dict.lookup("pistonCellSetName")),
    bowlInPistonPatchID_
    (
        dict.lookup("bowlInPistonPatchName"),
        mesh.boundaryMesh()
    ),
    bowlInCylinderPatchID_
    (
        dict.lookup("bowlInCylinderPatchName"),
        mesh.boundaryMesh()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::enginePiston::writeDict(Ostream& os) const
{
    simpleEnginePiston::writeDict(os);
    os  << nl << token::BEGIN_BLOCK
        << "pistonPointSetName " << pistonPointSetName_
        << token::END_STATEMENT << nl
        << "pistonFaceSetName " << pistonFaceSetName_
        << token::END_STATEMENT << nl
        << "pistonCellSetName " << pistonCellSetName_
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
