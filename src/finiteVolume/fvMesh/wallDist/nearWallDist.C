/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "nearWallDist.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "cellDistFuncs.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDist::doAll()
{
    const fvPatchList& patches = mesh_.boundary();

    forAll (patches, patchI)
    {
        fvPatchScalarField& yPatch = operator[](patchI);

        if (patches[patchI].isWall())
        {
            yPatch = 1/patches[patchI].deltaCoeffs();
        }
        else
        {
            yPatch = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    volScalarField::GeometricBoundaryField
    (
        mesh.boundary(),
        mesh.V(),           // Dummy internal field,
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(mesh)
{
    doAll();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallDist::correct()
{
    if (mesh_.changing())
    {
        // Update size if not equal
        if (size() != mesh_.boundary().size())
        {
            setSize(mesh_.boundary().size());
        }

        // Update size of GeometricBoundaryField
        forAll (mesh_.boundary(), patchI)
        {
            if (!set(patchI))
            {
                set
                (
                    patchI,
                    new calculatedFvPatchScalarField
                    (
                        mesh_.boundary()[patchI],
                        mesh_.V()                // Dummy internal field
                    )
                );
            }

            operator[](patchI).setSize(mesh_.boundary()[patchI].size());
        }
    }

    doAll();
}


// ************************************************************************* //
