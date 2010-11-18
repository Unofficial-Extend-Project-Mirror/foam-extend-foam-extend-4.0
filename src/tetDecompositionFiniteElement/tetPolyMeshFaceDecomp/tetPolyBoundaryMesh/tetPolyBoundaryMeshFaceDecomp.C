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

#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "polyBoundaryMesh.H"
#include "faceTetPolyPatchFaceDecomp.H"
#include "globalTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyBoundaryMesh
tetPolyBoundaryMeshFaceDecomp::tetPolyBoundaryMeshFaceDecomp
(
    const tetPolyMeshFaceDecomp& m,
    const polyBoundaryMesh& basicBdry
)
:
    tetPolyPatchFaceDecompList(basicBdry.size()),
    mesh_(m)
{
    // Set boundary patches
    tetPolyPatchFaceDecompList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.set
        (
            patchI,
            faceTetPolyPatchFaceDecomp::New(basicBdry[patchI], *this).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduInterfacePtrsList tetPolyBoundaryMeshFaceDecomp::interfaces() const
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


const globalTetPolyPatchFaceDecomp&
tetPolyBoundaryMeshFaceDecomp::globalPatch() const
{
    const tetPolyPatchFaceDecompList& patches = *this;

    forAll (patches, patchI)
    {
        if (isA<globalTetPolyPatchFaceDecomp>(patches[patchI]))
        {
            return refCast<const globalTetPolyPatchFaceDecomp>
            (
                patches[patchI]
            );
        }
    }

    FatalErrorIn
    (
        "const globalTetPolyPatchFaceDecomp&"
        "tetPolyBoundaryMeshFaceDecomp::globalPatch() const"
    )   << "patch not found.  Is this case running in parallel?"
        << abort(FatalError);

    // Dummy return
    return refCast<const globalTetPolyPatchFaceDecomp>(patches[0]);
}


faceListList tetPolyBoundaryMeshFaceDecomp::boundaryTriFaces() const
{
    faceListList result(size());

    forAll (result, patchI)
    {
        result[patchI] = operator[](patchI).triFaces();
    }

    return result;
}


void tetPolyBoundaryMeshFaceDecomp::updateMesh()
{
    tetPolyPatchFaceDecompList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches[patchI].updateMesh();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
