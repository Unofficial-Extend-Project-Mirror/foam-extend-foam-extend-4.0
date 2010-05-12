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

\*---------------------------------------------------------------------------*/

#include "faceTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<faceTetPolyPatchFaceDecomp> faceTetPolyPatchFaceDecomp::New
(
    const polyPatch& patch,
    const tetPolyBoundaryMeshFaceDecomp& bm
)
{
    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::New(const polyPatch&, "
            << " const tetPolyBoundaryMeshFaceDecomp&) : "
            << "constructing faceTetPolyPatchFaceDecomp"
            << endl;
    }

    polyPatchConstructorTable::iterator cstrIter =
        polyPatchConstructorTablePtr_->find(patch.type());

    if (cstrIter == polyPatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "faceTetPolyPatchFaceDecomp::New(const polyPatch&, "
            "const tetPolyBoundaryMeshFaceDecomp&)"
        )   << "Unknown faceTetPolyPatchFaceDecomp type "
            << patch.type()
            << ".  Valid faceTetPolyPatchFaceDecomp types are :" << endl
            << polyPatchConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<faceTetPolyPatchFaceDecomp>(cstrIter()(patch, bm));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
