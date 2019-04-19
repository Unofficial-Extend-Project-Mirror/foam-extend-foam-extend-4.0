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

#include "processorFaMeshes.H"
#include "foamTime.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorFaMeshes::readMeshes()
{
    forAll (fvMeshes_, procI)
    {
        if (fvMeshes_.set(procI))
        {
            meshes_.set
            (
                procI,
                new faMesh(fvMeshes_[procI])
            );

            pointProcAddressing_.set
           (
               procI,
               new labelIOList
               (
                   IOobject
                   (
                       "pointProcAddressing",
                       meshes_[procI].time().findInstance
                       (
                           meshes_[procI].meshDir(),
                           "pointProcAddressing"
                       ),
                       meshes_[procI].meshSubDir,
                       fvMeshes_[procI],
                       IOobject::MUST_READ,
                       IOobject::NO_WRITE
                   )
               )
           );

            edgeProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "edgeProcAddressing",
                        meshes_[procI].time().findInstance
                        (
                            meshes_[procI].meshDir(),
                            "edgeProcAddressing"
                        ),
                        meshes_[procI].meshSubDir,
                        fvMeshes_[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );

            faceProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "faceProcAddressing",
                        meshes_[procI].time().findInstance
                        (
                            meshes_[procI].meshDir(),
                            "faceProcAddressing"
                        ),
                        meshes_[procI].meshSubDir,
                        fvMeshes_[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );

            boundaryProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "boundaryProcAddressing",
                        meshes_[procI].time().findInstance
                        (
                            meshes_[procI].meshDir(),
                            "faceProcAddressing"
                        ),
                        meshes_[procI].meshSubDir,
                        fvMeshes_[procI],
                        IOobject::MUST_READ,
                    IOobject::NO_WRITE
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorFaMeshes::processorFaMeshes
(
    const PtrList<fvMesh>& processorFvMeshes,
    const bool read
)
:
    fvMeshes_(processorFvMeshes),
    meshes_(processorFvMeshes.size()),
    pointProcAddressing_(processorFvMeshes.size()),
    edgeProcAddressing_(processorFvMeshes.size()),
    faceProcAddressing_(processorFvMeshes.size()),
    boundaryProcAddressing_(processorFvMeshes.size())
{
    if (read)
    {
        readMeshes();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
