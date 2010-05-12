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

#include "processorFaMeshes.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorFaMeshes::read()
{
    forAll (fvMeshes_, procI)
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorFaMeshes::processorFaMeshes
(
    const PtrList<fvMesh>& processorFvMeshes
)
:
    fvMeshes_(processorFvMeshes),
    meshes_(processorFvMeshes.size()),
    pointProcAddressing_(processorFvMeshes.size()),
    edgeProcAddressing_(processorFvMeshes.size()),
    faceProcAddressing_(processorFvMeshes.size()),
    boundaryProcAddressing_(processorFvMeshes.size())
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Foam::fvMesh::readUpdateState Foam::processorFaMeshes::readUpdate()
// {
//     fvMesh::readUpdateState stat = fvMesh::UNCHANGED;

//     forAll (databases_, procI)
//     {
//         // Check if any new meshes need to be read.
//         fvMesh::readUpdateState procStat = meshes_[procI].readUpdate();

//         /*
//         if (procStat != fvMesh::UNCHANGED)
//         {
//             Info<< "Processor " << procI
//                 << " at time " << databases_[procI].timeName()
//                 << " detected mesh change " << procStat
//                 << endl;
//         }
//         */

//         // Combine into overall mesh change status
//         if (stat == fvMesh::UNCHANGED)
//         {
//             stat = procStat;
//         }
//         else
//         {
//             if (stat != procStat)
//             {
//                 FatalErrorIn("processorFaMeshes::readUpdate()")
//                     << "Processor " << procI
//                     << " has a different polyMesh at time "
//                     << databases_[procI].timeName()
//                     << " compared to any previous processors." << nl
//                     << "Please check time " << databases_[procI].timeName()
//                     << " directories on all processors for consistent"
//                     << " mesh files."
//                     << exit(FatalError);
//             }
//         }
//     }

//     if
//     (
//         stat == fvMesh::TOPO_CHANGE
//      || stat == fvMesh::TOPO_PATCH_CHANGE
//     )
//     {
//         // Reread all meshes and addresssing
//         read();
//     }
//     return stat;
// }


// void Foam::processorFaMeshes::reconstructPoints(fvMesh& mesh)
// {
//     // Read the field for all the processors
//     PtrList<pointIOField> procsPoints(meshes_.size());

//     forAll (meshes_, procI)
//     {
//         procsPoints.set
//         (
//             procI,
//             new pointIOField
//             (
//                 IOobject
//                 (
//                     "points",
//                     meshes_[procI].time().timeName(),
//                     polyMesh::meshSubDir,
//                     meshes_[procI],
//                     IOobject::MUST_READ,
//                     IOobject::NO_WRITE
//                 )
//             )
//         );
//     }

//     // Create the new points
//     vectorField newPoints(mesh.nPoints());

//     forAll (meshes_, procI)
//     {
//         const vectorField& procPoints = procsPoints[procI];

//         // Set the cell values in the reconstructed field

//         const labelList& pointProcAddressingI = pointProcAddressing_[procI];

//         if (pointProcAddressingI.size() != procPoints.size())
//         {
//             FatalErrorIn("processorFaMeshes")
//                 << "problem :"
//                 << " pointProcAddressingI:" << pointProcAddressingI.size()
//                 << " procPoints:" << procPoints.size()
//                 << abort(FatalError);
//         }

//         forAll(pointProcAddressingI, pointI)
//         {
//             newPoints[pointProcAddressingI[pointI]] = procPoints[pointI];
//         }
//     }

//     mesh.movePoints(newPoints);
//     mesh.write();
// }


// ************************************************************************* //
