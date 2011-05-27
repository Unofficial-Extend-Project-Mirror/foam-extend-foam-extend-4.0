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

#include "processorMeshesReconstructor.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorMeshesReconstructor::readMeshes()
{
    forAll (databases_, procI)
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
                    databases_[procI].timeName(),
                    databases_[procI]
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
    PtrList<Time>& databases,
    const word& meshName
)
:
    databases_(databases),
    meshName_(meshName),
    meshes_(databases.size()),
    pointProcAddressing_(),
    faceProcAddressing_(),
    cellProcAddressing_(),
    boundaryProcAddressing_()
{
    readMeshes();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState
Foam::processorMeshesReconstructor::readUpdate()
{
    polyMesh::readUpdateState stat = polyMesh::UNCHANGED;

    forAll (databases_, procI)
    {
        // Check if any new meshes need to be read.
        polyMesh::readUpdateState procStat = meshes_[procI].readUpdate();

        // Combine into overall mesh change status
        if (stat == polyMesh::UNCHANGED)
        {
            stat = procStat;
        }
        else
        {
            if (stat != procStat)
            {
                FatalErrorIn("processorMeshesReconstructor::readUpdate()")
                    << "Processor " << procI
                    << " has a different polyMesh at time "
                    << databases_[procI].timeName()
                    << " compared to any previous processors." << nl
                    << "Please check time " << databases_[procI].timeName()
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
