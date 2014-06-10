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

Description

\*---------------------------------------------------------------------------*/

#include "polyMeshGen.H"
#include "demandDrivenData.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGen::polyMeshGen(const Time& t)
:
    polyMeshGenCells(t),
    metaDict_()
{}

//- Construct from components without the boundary
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    polyMeshGenCells(t, points, faces, cells),
    metaDict_()
{}

//- Construct from components with the boundary
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const wordList& patchNames,
    const labelList& patchStart,
    const labelList& nFacesInPatch
)
:
    polyMeshGenCells
    (
        t,
        points,
        faces,
        cells,
        patchNames,
        patchStart,
        nFacesInPatch
    ),
    metaDict_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
polyMeshGen::~polyMeshGen()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGen::read()
{
    polyMeshGenCells::read();

    metaDict_ =
        IOdictionary
        (
            IOobject
            (
                "meshMetaDict",
                runTime_.constant(),
                "polyMesh",
                runTime_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );
}

void polyMeshGen::write() const
{
    //- remove old mesh before writting
    const fileName meshDir = runTime_.path()/runTime_.constant()/"polyMesh";

    rm(meshDir/"points");
    rm(meshDir/"faces");
    rm(meshDir/"owner");
    rm(meshDir/"neighbour");
    rm(meshDir/"cells");
    rm(meshDir/"boundary");
    rm(meshDir/"pointZones");
    rm(meshDir/"faceZones");
    rm(meshDir/"cellZones");
    rm(meshDir/"meshModifiers");
    rm(meshDir/"parallelData");
    rm(meshDir/"meshMetaDict");

    // remove sets if they exist
    if (isDir(meshDir/"sets"))
    {
        rmDir(meshDir/"sets");
    }

    //- write the mesh
    polyMeshGenCells::write();

    //- write meta data
    OFstream fName(meshDir/"meshMetaDict");
    IOdictionary writeMeta
    (
        IOobject
        (
            "meshMetaDict",
            runTime_.constant(),
            "polyMesh",
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        metaDict_
    );

    writeMeta.writeHeader(fName);
    writeMeta.writeData(fName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
