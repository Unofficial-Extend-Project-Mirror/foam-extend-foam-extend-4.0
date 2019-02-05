/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
    edgesPtr_(nullptr),
    ccPtr_(nullptr),
    ecPtr_(nullptr),
    pcPtr_(nullptr),
    efPtr_(nullptr),
    pfPtr_(nullptr),

    cePtr_(nullptr),
    fePtr_(nullptr),
    pePtr_(nullptr),
    ppPtr_(nullptr),
    cpPtr_(nullptr),
    cellCentresPtr_(nullptr),
    faceCentresPtr_(nullptr),
    cellVolumesPtr_(nullptr),
    faceAreasPtr_(nullptr),
    globalPointLabelPtr_(nullptr),
    globalFaceLabelPtr_(nullptr),
    globalCellLabelPtr_(nullptr),
    globalEdgeLabelPtr_(nullptr),
    pProcsPtr_(nullptr),
    globalToLocalPointAddressingPtr_(nullptr),
    pointNeiProcsPtr_(nullptr),
    eProcsPtr_(nullptr),
    globalToLocalEdgeAddressingPtr_(nullptr),
    edgeNeiProcsPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMeshGenAddressing::~polyMeshGenAddressing()
{
    clearAll();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
