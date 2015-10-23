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

Class
    faGlobalMeshData

Description

Author
    Hrvoje Jasak

\*----------------------------------------------------------------------------*/

#include "faGlobalMeshData.H"
#include "faMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faGlobalMeshData::faGlobalMeshData(const faMesh& mesh)
:
    faProcessorTopology(mesh.boundary()),
    mesh_(mesh),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faGlobalMeshData::~faGlobalMeshData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faMesh& Foam::faGlobalMeshData::mesh() const
{
    return mesh_;
}


// Update all data after morph
void Foam::faGlobalMeshData::updateMesh()
{
    label polyMeshNGlobalPoints =
        mesh_().globalData().nGlobalPoints();

    const labelList& polyMeshSharedPointLabels =
        mesh_().globalData().sharedPointLabels();

    const labelList& polyMeshSharedPointAddr =
        mesh_().globalData().sharedPointAddr();

    labelHashSet sharedPointLabels;

    labelField globalList(polyMeshNGlobalPoints, 0);

    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].type() == processorFaPatch::typeName)
        {
            const labelList& localPointLabels =
                mesh_.boundary()[patchI].pointLabels();

            forAll(localPointLabels, pointI)
            {
                label polyMeshPoint =
                    mesh_.patch().meshPoints()[localPointLabels[pointI]];

                label sharedPolyMeshPoint =
                    findIndex(polyMeshSharedPointLabels, polyMeshPoint);

                if
                (
                    sharedPolyMeshPoint != -1
                 && !sharedPointLabels.found(localPointLabels[pointI])
                )
                {
                    globalList[polyMeshSharedPointAddr[sharedPolyMeshPoint]]
                        += 1;

                    sharedPointLabels.insert(localPointLabels[pointI]);
                }
            }
        }
    }

    sharedPointLabels_ = sharedPointLabels.toc();

    combineReduce(globalList, plusEqOp<labelField >());

    nGlobalPoints_ = 0;
    for (label i=0; i<globalList.size(); i++)
    {
        if(globalList[i] > 0)
        {
            globalList[i] = ++nGlobalPoints_;
        }
    }

    sharedPointAddr_.setSize(sharedPointLabels_.size());
    forAll(sharedPointAddr_, pointI)
    {
        label polyMeshSharedPointIndex = findIndex
        (
            polyMeshSharedPointLabels,
            mesh_.patch().meshPoints()[sharedPointLabels_[pointI]]
        );

        sharedPointAddr_[pointI] =
            globalList[polyMeshSharedPointAddr[polyMeshSharedPointIndex]]
          - 1;
    }
}

// ************************************************************************* //
