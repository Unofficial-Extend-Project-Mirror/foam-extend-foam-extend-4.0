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

#include "processorMeshesReconstructor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::processorMeshesReconstructor, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorMeshesReconstructor::readMeshes(PtrList<Time>& databases)
{
    forAll (databases, procI)
    {
        Info<< "Reading mesh for processor " << procI << endl;
        meshes_.set
        (
            procI,
            new fvMesh
            (
                IOobject
                (
                    meshName_,
                    databases[procI].timeName(),
                    databases[procI]
                )
            )
        );
    }

    // Clear reconstruction maps: new mesh
    clearMaps();
}


void Foam::processorMeshesReconstructor::clearMaps()
{
    pointProcAddressing_.clear();
    faceProcAddressing_.clear();
    cellProcAddressing_.clear();
    boundaryProcAddressing_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::processorMeshesReconstructor::processorMeshesReconstructor
(
    const word& meshName
)
:
    meshName_(meshName),
    meshes_(),
    globalPointIndex_(),
    pointProcAddressing_(),
    faceProcAddressing_(),
    cellProcAddressing_(),
    boundaryProcAddressing_()
{}


Foam::processorMeshesReconstructor::processorMeshesReconstructor
(
    PtrList<Time>& databases,
    const word& meshName
)
:
    meshName_(meshName),
    meshes_(databases.size()),
    globalPointIndex_(),
    pointProcAddressing_(),
    faceProcAddressing_(),
    cellProcAddressing_(),
    boundaryProcAddressing_()
{
    readMeshes(databases);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState
Foam::processorMeshesReconstructor::readUpdate()
{
    polyMesh::readUpdateState stat = polyMesh::UNCHANGED;

    forAll (meshes_, procI)
    {
        // Only do action if database has been set
        if (meshes_.set(procI))
        {
            // Check if any new meshes need to be read.
            polyMesh::readUpdateState procStat = meshes_[procI].readUpdate();

            // Combine into overall mesh change status
            if (stat == polyMesh::UNCHANGED)
            {
                stat = procStat;
            }
            else if (stat != procStat)
            {
                if (
                    (
                        stat == polyMesh::TOPO_CHANGE
                     && procStat == polyMesh::TOPO_PATCH_CHANGE
                    )
                 || (
                        procStat == polyMesh::TOPO_CHANGE
                     && stat == polyMesh::TOPO_PATCH_CHANGE
                    )
                )
                {
                    continue;
                }

                FatalErrorIn("processorMeshesReconstructor::readUpdate()")
                    << "Processor " << procI
                    << " has a different polyMesh at time "
                    << meshes_[procI].time().timeName()
                    << " compared to any previous processors." << nl
                    << "Please check time "
                    << meshes_[procI].time().timeName()
                    << " directories on all processors for consistent"
                    << " mesh files."
                    << exit(FatalError);
            }
        }
    }

    if
    (
        stat == polyMesh::TOPO_CHANGE
     || stat == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        // Reconstruction maps have changed: clear
        Info<< "Mesh changed: clear reconstruction maps" << endl;
        clearMaps();
    }

    return stat;
}


const Foam::PtrList<Foam::labelIOList>&
Foam::processorMeshesReconstructor::pointProcAddressing() const
{
    if (pointProcAddressing_.empty())
    {
        FatalErrorIn
        (
            "const PtrList<labelIOList>& "
            "processorMeshesReconstructor::pointProcAddressing() const"
        )   << "Mesh is not reconstructed"
            << abort(FatalError);
    }

    return pointProcAddressing_;
}


const Foam::PtrList<Foam::labelIOList>&
Foam::processorMeshesReconstructor::faceProcAddressing() const
{
    if (faceProcAddressing_.empty())
    {
        FatalErrorIn
        (
            "const PtrList<labelIOList>& "
            "processorMeshesReconstructor::faceProcAddressing() const"
        )   << "Mesh is not reconstructed"
            << abort(FatalError);
    }

    return faceProcAddressing_;
}


const Foam::PtrList<Foam::labelIOList>&
Foam::processorMeshesReconstructor::cellProcAddressing() const
{
    if (cellProcAddressing_.empty())
    {
        FatalErrorIn
        (
            "const PtrList<labelIOList>& "
            "processorMeshesReconstructor::cellProcAddressing() const"
        )   << "Mesh is not reconstructed"
            << abort(FatalError);
    }

    return cellProcAddressing_;
}


const Foam::PtrList<Foam::labelIOList>&
Foam::processorMeshesReconstructor::boundaryProcAddressing() const
{
    if (boundaryProcAddressing_.empty())
    {
        FatalErrorIn
        (
            "const PtrList<labelIOList>& "
            "processorMeshesReconstructor::boundaryProcAddressing() const"
        )   << "Mesh is not reconstructed"
            << abort(FatalError);
    }

    return boundaryProcAddressing_;
}


// ************************************************************************* //
