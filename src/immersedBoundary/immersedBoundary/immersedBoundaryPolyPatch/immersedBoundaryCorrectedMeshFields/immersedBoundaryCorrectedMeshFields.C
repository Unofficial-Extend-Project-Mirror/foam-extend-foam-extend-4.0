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

#include "immersedBoundaryCorrectedMeshFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryCorrectedMeshFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::immersedBoundaryCorrectedMeshFields::immersedBoundaryCorrectedMeshFields
(
    const polyMesh& mesh
)
:
    MeshObject<polyMesh, immersedBoundaryCorrectedMeshFields>(mesh),
    correctedCellCentresPtr_(nullptr),
    correctedFaceCentresPtr_(nullptr),
    correctedCellVolumesPtr_(nullptr),
    correctedFaceAreasPtr_(nullptr),
    curTopoIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::immersedBoundaryCorrectedMeshFields::~immersedBoundaryCorrectedMeshFields()
{
    deleteDemandDrivenData(correctedCellCentresPtr_);
    deleteDemandDrivenData(correctedFaceCentresPtr_);
    deleteDemandDrivenData(correctedCellVolumesPtr_);
    deleteDemandDrivenData(correctedFaceAreasPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedBoundaryCorrectedMeshFields::clearOut
(
    const label topoChangeIndex
) const
{
    // Check if the mesh changed since the last call and delete data if so.
    // The fields must not get reset by every IB poly patch, this is why this
    // check is necessary

    if (curTopoIndex_ < topoChangeIndex)
    {
        deleteDemandDrivenData(correctedCellCentresPtr_);
        deleteDemandDrivenData(correctedFaceCentresPtr_);
        deleteDemandDrivenData(correctedCellVolumesPtr_);
        deleteDemandDrivenData(correctedFaceAreasPtr_);

        curTopoIndex_ = topoChangeIndex;
    }
}


Foam::vectorField&
Foam::immersedBoundaryCorrectedMeshFields::correctedCellCentres() const
{
    if (!correctedCellCentresPtr_)
    {
        correctedCellCentresPtr_ = new vectorField(mesh().cellCentres());
    }

    return *correctedCellCentresPtr_;
}


Foam::vectorField&
Foam::immersedBoundaryCorrectedMeshFields::correctedFaceCentres() const
{
    if (!correctedFaceCentresPtr_)
    {
        correctedFaceCentresPtr_ = new vectorField(mesh().faceCentres());
    }

    return *correctedFaceCentresPtr_;
}


Foam::scalarField&
Foam::immersedBoundaryCorrectedMeshFields::correctedCellVolumes() const
{
    if (!correctedCellVolumesPtr_)
    {
        correctedCellVolumesPtr_ = new scalarField(mesh().cellVolumes());
    }

    return *correctedCellVolumesPtr_;
}


Foam::vectorField&
Foam::immersedBoundaryCorrectedMeshFields::correctedFaceAreas() const
{
    if (!correctedFaceAreasPtr_)
    {
        correctedFaceAreasPtr_ = new vectorField(mesh().faceAreas());
    }

    return *correctedFaceAreasPtr_;
}


// ************************************************************************* //
