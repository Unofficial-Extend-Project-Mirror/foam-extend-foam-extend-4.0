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

Description

\*---------------------------------------------------------------------------*/

#include "vtkMesh.H"
#include "fvMeshSubset.H"
#include "foamTime.H"
#include "cellSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::vtkMesh::vtkMesh
(
    const IOobject& io,
    const word& setName
)
:
    fvMesh(io),
    subsetMesh_
    (
        IOobject
        (
            "subset",
            io.time().constant(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    setName_(setName)
{
    if (setName.size())
    {
        // Read cellSet using whole mesh
        cellSet currentSet(*this, setName_);

        // Set current subset
        subsetMesh_.setLargeCellSubset(currentSet);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkMesh::~vtkMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState Foam::vtkMesh::readUpdate()
{
    polyMesh::readUpdateState meshState = fvMesh::readUpdate();

    if (meshState != polyMesh::UNCHANGED)
    {
        // Note: since fvMeshSubset has no movePoints() functionality,
        //  reconstruct the subset even if only movement.

        topoPtr_.clear();

        if (setName_.size())
        {
            Info<< "Subsetting mesh based on cellSet " << setName_ << endl;

            // Read cellSet using whole mesh
            cellSet currentSet(*this, setName_);

            subsetMesh_.setLargeCellSubset(currentSet);
        }
    }

    return meshState;
}


// ************************************************************************* //
