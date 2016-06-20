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

Description

\*---------------------------------------------------------------------------*/

#include "MeshWave.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamReduceOps.H"
#include "debug.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Iterate, propagating changedFacesInfo across mesh, until no change (or
// maxIter reached).
template <class Type>
Foam::MeshWave<Type>::MeshWave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    const label maxIter
)
:
    allFaceInfo_(mesh.nFaces()),
    allCellInfo_(mesh.nCells()),
    calc_
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        allFaceInfo_,
        allCellInfo_,
        maxIter
    )
{}


// Iterate, propagating changedFacesInfo across mesh, until no change (or
// maxIter reached). Initial cell values specified.
template <class Type>
Foam::MeshWave<Type>::MeshWave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    const List<Type>& allCellInfo,
    const label maxIter
)
:
    allFaceInfo_(mesh.nFaces()),
    allCellInfo_(allCellInfo),
    calc_
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        allFaceInfo_,
        allCellInfo_,
        maxIter
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
