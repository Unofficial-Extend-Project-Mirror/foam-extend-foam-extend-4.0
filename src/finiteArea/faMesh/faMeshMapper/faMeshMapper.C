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

#include "faMeshMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshMapper::faMeshMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    nOldPoints_(mesh.nPoints()),
    nOldEdges_(mesh.nEdges()),
    nOldInternalEdges_(mesh.nInternalEdges()),
    nOldFaces_(mesh.nFaces()),
    oldPatchSizes_(mesh.boundary().size(), 0),
    oldPatchStarts_(mesh.boundary().size(), -1),
    oldPatchEdgeFaces_(mesh.boundary().size()),
    areaMap_(mesh, mpm),
    edgeMap_(mesh, mpm),
    boundaryMap_(mesh, mpm)
{
    // Capture old patch information
    const faBoundaryMesh& patches = mesh.boundary();

    forAll (patches, patchI)
    {
        oldPatchSizes_[patchI] = patches[patchI].size();
        oldPatchStarts_[patchI] = patches[patchI].start();

        oldPatchEdgeFaces_[patchI] = patches[patchI].edgeFaces();
    }
}


// ************************************************************************* //
