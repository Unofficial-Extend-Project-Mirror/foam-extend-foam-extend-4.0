/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "reflectionVectors.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::reflectionVectors::reflectionVectors(const Foam::fvMesh& mesh)
:
    n_
    (
        IOobject
        (
            "reflectionVectors",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("n", dimless, vector::zero)
    )
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reflectionVectors::correct()
{
    const fvMesh& mesh = n_.mesh();
    const fvPatchList& patches = mesh.boundary();

    forAll (patches, patchi)
    {
        // find the nearest face for every cell
        if (patches[patchi].isWall())
        {
            n_.boundaryField()[patchi] =
                mesh.Sf().boundaryField()[patchi]
               /mesh.magSf().boundaryField()[patchi];
        }
    }
}


// ************************************************************************* //
