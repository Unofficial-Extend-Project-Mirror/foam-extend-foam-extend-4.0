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

#include "faceTetPolyPatchCellDecomp.H"
#include "tetPolyBoundaryMeshCellDecomp.H"
#include "tetPolyMeshCellDecomp.H"
#include "demandDrivenData.H"
#include "boolList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faceTetPolyPatchCellDecomp, 0);
defineRunTimeSelectionTable(faceTetPolyPatchCellDecomp, polyPatch);

addToRunTimeSelectionTable
(
    faceTetPolyPatchCellDecomp,
    faceTetPolyPatchCellDecomp,
    polyPatch
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyPatch
faceTetPolyPatchCellDecomp::faceTetPolyPatchCellDecomp
(
    const polyPatch& p,
    const tetPolyBoundaryMeshCellDecomp& bm
)
:
    tetPolyPatchCellDecomp(bm),
    boundaryIndex_(p.index())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const polyPatch& faceTetPolyPatchCellDecomp::patch() const
{
    return boundaryMesh().mesh()().polyMesh::boundaryMesh()[index()];
}


triFaceList faceTetPolyPatchCellDecomp::faceTriangles
(
    const label faceID
) const
{
    const face& f =
        boundaryMesh().mesh()()
            .boundaryMesh()[index()].localFaces()[faceID];

    // Create a list of triangles to keep the triangles that
    // have already been added
    triFaceList result(f.size() - 2);

    for (label triI = 0; triI < (f.size() - 2); triI++)
    {
        result[triI] = triFace(f[0], f[triI + 1], f[triI + 2]);
    }

    return result;
}


faceList faceTetPolyPatchCellDecomp::triFaces() const
{
    const faceList& mf =
        boundaryMesh().mesh()().boundaryMesh()[index()];

    // Count the number of faces
    label nTriFaces = 0;

    forAll (mf, faceI)
    {
        nTriFaces += mf[faceI].size() - 2;
    }

    faceList result(nTriFaces);

    // Reset the counter for re-use
    nTriFaces = 0;

    face triangleFace(3);

    forAll (mf, faceI)
    {
        const face& f = mf[faceI];

        for (label triI = 0; triI < (f.size() - 2); triI++)
        {
            triangleFace[0] = f[0];
            triangleFace[1] = f[triI + 1];
            triangleFace[2] = f[triI + 2];

            result[nTriFaces] = triangleFace;
            nTriFaces++;
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
