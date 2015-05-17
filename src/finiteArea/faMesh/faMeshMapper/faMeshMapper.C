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
