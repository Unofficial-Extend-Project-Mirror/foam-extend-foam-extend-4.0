/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyMeshGenAddressing, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyMeshGenAddressing::polyMeshGenAddressing(const polyMeshGenCells& mesh)
:
    mesh_(mesh),
    edgesPtr_(NULL),
    ccPtr_(NULL),
    ecPtr_(NULL),
    pcPtr_(NULL),
    efPtr_(NULL),
    pfPtr_(NULL),

    cePtr_(NULL),
    fePtr_(NULL),
    pePtr_(NULL),
    ppPtr_(NULL),
    cpPtr_(NULL),
    cellCentresPtr_(NULL),
    faceCentresPtr_(NULL),
    cellVolumesPtr_(NULL),
    faceAreasPtr_(NULL),
    globalPointLabelPtr_(NULL),
    globalFaceLabelPtr_(NULL),
    globalCellLabelPtr_(NULL),
    globalEdgeLabelPtr_(NULL),
    pProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    pointNeiProcsPtr_(NULL),
    eProcsPtr_(NULL),
    globalToLocalEdgeAddressingPtr_(NULL),
    edgeNeiProcsPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMeshGenAddressing::~polyMeshGenAddressing()
{
    clearAll();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
