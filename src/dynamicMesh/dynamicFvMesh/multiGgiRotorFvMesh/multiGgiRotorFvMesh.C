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

#include "multiGgiRotorFvMesh.H"
#include "foamTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiGgiRotorFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, multiGgiRotorFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::multiGgiRotorFvMesh::multiGgiRotorFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dict_
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
    ),
    rotors_()
{
    // Read rotors from the dictionary
    PtrList<entry> rotorEntries(dict_.lookup("rotors"));
    rotors_.setSize(rotorEntries.size());

    forAll (rotorEntries, rotorI)
    {
        rotors_.set
        (
            rotorI,
            new ggiRotor
            (
                rotorEntries[rotorI].keyword(),
                *this,
                rotorEntries[rotorI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiGgiRotorFvMesh::~multiGgiRotorFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiGgiRotorFvMesh::update()
{
    forAll (rotors_, rotorI)
    {
        rotors_[rotorI].updateTopology();
    }

    // Accumulate point motion
    vectorField pointMotion(allPoints().size(), vector::zero);

    forAll (rotors_, rotorI)
    {
        pointMotion += rotors_[rotorI].pointMotion();
    }

    // Move points
    movePoints(allPoints() + pointMotion);

    // Motion only - return false
    return true;
}


// ************************************************************************* //
