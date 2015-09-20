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

#include "pointFieldDecomposer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const pointPatch& completeMeshPatch,
    const pointPatch& procMeshPatch,
    const labelList& directAddr
)
:
    PointPatchFieldMapperPatchRef<pointPatch>
    (
        completeMeshPatch,
        procMeshPatch
    ),
    sizeBeforeMapping_(completeMeshPatch.size()),
    directAddressing_(procMeshPatch.size(), -1)
{
    // Create the inverse-addressing of the patch point labels.
    labelList pointMap(completeMeshPatch.boundaryMesh().mesh().size(), -1);

    const labelList& completeMeshPatchPoints = completeMeshPatch.meshPoints();

    forAll (completeMeshPatchPoints, pointi)
    {
        pointMap[completeMeshPatchPoints[pointi]] = pointi;
    }

    // Use the inverse point addressing to create the addressing table for this
    // patch
    const labelList& procMeshPatchPoints = procMeshPatch.meshPoints();

    forAll (procMeshPatchPoints, pointi)
    {
        directAddressing_[pointi] =
            pointMap[directAddr[procMeshPatchPoints[pointi]]];
    }

    // Check that all the patch point addresses are set
    if (directAddressing_.size() && min(directAddressing_) < 0)
    {
        FatalErrorIn
        (
            "pointFieldDecomposer::patchFieldDecomposer()"
        )   << "Incomplete patch point addressing"
            << abort(FatalError);
    }
}


pointFieldDecomposer::pointFieldDecomposer
(
    const pointMesh& completeMesh,
    const pointMesh& procMesh,
    const labelList& pointAddressing,
    const labelList& boundaryAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    pointAddressing_(pointAddressing),
    boundaryAddressing_(boundaryAddressing),
    patchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<patchFieldDecomposer*>(NULL)
    )
{
    forAll (boundaryAddressing_, patchi)
    {
        if (boundaryAddressing_[patchi] >= 0)
        {
            patchFieldDecomposerPtrs_[patchi] = new patchFieldDecomposer
            (
                completeMesh_.boundary()[boundaryAddressing_[patchi]],
                procMesh_.boundary()[patchi],
                pointAddressing_
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pointFieldDecomposer::~pointFieldDecomposer()
{
    forAll (patchFieldDecomposerPtrs_, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            delete patchFieldDecomposerPtrs_[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
