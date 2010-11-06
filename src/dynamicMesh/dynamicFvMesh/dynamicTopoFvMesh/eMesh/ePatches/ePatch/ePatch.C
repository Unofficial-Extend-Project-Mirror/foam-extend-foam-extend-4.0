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

Description

\*---------------------------------------------------------------------------*/

#include "ePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "eBoundaryMesh.H"
#include "eMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ePatch, 0);

defineRunTimeSelectionTable(ePatch, word);
defineRunTimeSelectionTable(ePatch, dictionary);

addToRunTimeSelectionTable(ePatch, ePatch, word);
addToRunTimeSelectionTable(ePatch, ePatch, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ePatch::ePatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const eBoundaryMesh& bm
)
:
    patchIdentifier(name, index),
    boundaryMesh_(bm),
    start_(start),
    size_(size)
{}


// Construct from dictionary
ePatch::ePatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const eBoundaryMesh& bm
)
:
    patchIdentifier(name, dict, index),
    boundaryMesh_(bm),
    start_(readLabel(dict.lookup("start"))),
    size_(readLabel(dict.lookup("size")))
{}


ePatch::ePatch(const ePatch& p)
:
    patchIdentifier(p, p.index()),
    boundaryMesh_(p.boundaryMesh_),
    start_(p.start()),
    size_(p.size())
{}


ePatch::ePatch(const ePatch& p, const eBoundaryMesh& bm)
:
    patchIdentifier(p, p.index()),
    boundaryMesh_(bm),
    start_(p.start()),
    size_(p.size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ePatch::~ePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const eBoundaryMesh& ePatch::boundaryMesh() const
{
    return boundaryMesh_;
}


void ePatch::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    os.writeKeyword("start") << start() << token::END_STATEMENT << nl;
    os.writeKeyword("size") << size() << token::END_STATEMENT << nl;
    patchIdentifier::write(os);
}


//- Assignment operator
void ePatch::operator=(const ePatch& p)
{
    patchIdentifier::operator=(p);
    start_ = p.start_;
    size_ = p.size_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const ePatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const faPatch& p)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
