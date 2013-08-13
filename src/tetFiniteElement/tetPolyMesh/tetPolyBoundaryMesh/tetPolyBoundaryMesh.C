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

#include "tetPolyBoundaryMesh.H"
#include "polyBoundaryMesh.H"
#include "faceTetPolyPatch.H"
#include "globalTetPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyBoundaryMesh
tetPolyBoundaryMesh::tetPolyBoundaryMesh
(
    const tetPolyMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    tetPolyPatchList(basicBdry.size()),
    mesh_(m)
{
    // Set boundary patches
    tetPolyPatchList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.set
        (
            patchI,
            faceTetPolyPatch::New(basicBdry[patchI], *this).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduInterfacePtrsList tetPolyBoundaryMesh::interfaces() const
{
    lduInterfacePtrsList interfaces(size());

    forAll (interfaces, patchi)
    {
        if (isA<lduInterface>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
               &refCast<const lduInterface>(this->operator[](patchi))
            );
        }
    }

    return interfaces;
}


const globalTetPolyPatch&
tetPolyBoundaryMesh::globalPatch() const
{
    const tetPolyPatchList& patches = *this;

    forAll (patches, patchI)
    {
        if (isA<globalTetPolyPatch>(patches[patchI]))
        {
            return refCast<const globalTetPolyPatch>
            (
                patches[patchI]
            );
        }
    }

    FatalErrorIn
    (
        "const globalTetPolyPatch&"
        "tetPolyBoundaryMesh::globalPatch() const"
    )   << "patch not found.  Is this case running in parallel?"
        << abort(FatalError);

    // Dummy return
    return refCast<const globalTetPolyPatch>(patches[0]);
}


faceListList tetPolyBoundaryMesh::boundaryTriFaces() const
{
    faceListList result(size());

    forAll (result, patchI)
    {
        result[patchI] = operator[](patchI).triFaces();
    }

    return result;
}


void tetPolyBoundaryMesh::updateMesh()
{
    tetPolyPatchList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches[patchI].updateMesh();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
