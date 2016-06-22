/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "attachDetachFvMesh.H"
#include "foamTime.H"
#include "attachDetach.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(attachDetachFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        attachDetachFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::attachDetachFvMesh::addZonesAndModifiers()
{
    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add a topology modifier
    Info << "Adding topology modifiers" << endl;

    PtrList<entry> entries(motionDict_.lookup("modifiers"));

    topoChanger_.setSize(entries.size());

    forAll (entries, i)
    {
        topoChanger_.set
        (
            i,
            new attachDetach
            (
                entries[i].keyword(),
                i,
                topoChanger_,
                word(entries[i].dict().lookup("faceZone")),
                word(entries[i].dict().lookup("topPatch")),
                word(entries[i].dict().lookup("bottomPatch")),
                scalarField(entries[i].dict().lookup("triggerTimes"))
            )
        );
    }

    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::attachDetachFvMesh::attachDetachFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    )
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::attachDetachFvMesh::~attachDetachFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::attachDetachFvMesh::update()
{
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    return topoChangeMap->morphing();
}


// ************************************************************************* //
