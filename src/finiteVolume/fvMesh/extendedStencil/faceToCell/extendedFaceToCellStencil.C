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

#include "extendedFaceToCellStencil.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedFaceToCellStencil::extendedFaceToCellStencil(const polyMesh& mesh)
:
    mesh_(mesh)
{
    // Check for transformation - not supported.
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<coupledPolyPatch>(patches[patchI]))
        {
            const coupledPolyPatch& cpp =
                refCast<const coupledPolyPatch>(patches[patchI]);

            if (!cpp.parallel() || cpp.separated())
            {
                FatalErrorIn
                (
                    "extendedFaceToCellStencil::extendedFaceToCellStencil"
                    "(const polyMesh&)"
                )   << "Coupled patches with transformations not supported."
                    << endl
                    << "Problematic patch " << cpp.name() << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
