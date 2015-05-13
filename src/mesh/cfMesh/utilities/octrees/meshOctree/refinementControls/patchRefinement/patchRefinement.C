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

#include "patchRefinement.H"
#include "triSurf.H"
#include "Istream.H"
#include "demandDrivenData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchRefinement::patchRefinement()
:
    patchName_(),
    cellSize_(0.0)
{
}

patchRefinement::patchRefinement(const word& n, const scalar& s)
:
    patchName_(n),
    cellSize_(s)
{
}

patchRefinement::patchRefinement(Istream& is)
:
    patchName_(is),
    cellSize_(readScalar(is))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

patchRefinement::~patchRefinement()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word patchRefinement::patchName() const
{
    return patchName_;
}

label patchRefinement::patchInSurface(const triSurf& ts) const
{
    forAll(ts.patches(), patchI)
        if( ts.patches()[patchI].name() == patchName_ )
            return patchI;

    FatalErrorIn
    (
        "label patchRefinement::patchInSurface(const triSurf& ts) const"
    ) << "Patch " << patchName_
        << " does not exist in surface " << ts.patches()
        << exit(FatalError);

    return -1;
}

scalar patchRefinement::cellSize() const
{
    return cellSize_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void patchRefinement::operator=(const patchRefinement& pr)
{
    patchName_ = pr.patchName_;
    cellSize_ = pr.cellSize_;
}

Istream& operator>>(Istream& is, patchRefinement& pr)
{
    pr.patchName_ = word(is);
    pr.cellSize_ = readScalar(is);
    return is;
}

Ostream& operator<<(Ostream& os, const patchRefinement& pr)
{
    os << pr.patchName_ << " " << pr.cellSize_ << nl;
    return os;
}

bool operator==(const patchRefinement& p1, const patchRefinement& p2)
{
    if( p1.patchName() == p2.patchName() )
        return true;

    return false;
}

bool operator!=(const patchRefinement& p1, const patchRefinement& p2)
{
    return !(p1 == p2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
