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

#include "boundaryPatch.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundaryPatch, 0);
addToRunTimeSelectionTable(boundaryPatchBase, boundaryPatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

boundaryPatch::boundaryPatch
(
    const word& n,
    const word& t,
    const label nF,
    const label sF
)
:
    boundaryPatchBase(n, t, nF, sF)
{}

boundaryPatch::boundaryPatch(const word& name, const dictionary& dict)
:
    boundaryPatchBase(name, dict)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary boundaryPatch::dict() const
{
    dictionary dict;

    dict.add("type", type_);

    dict.add("nFaces", nFaces_);
    dict.add("startFace", startFace_);

    return dict;
}

void boundaryPatch::write(Ostream& os) const
{
    this->operator<<(os);
}

void boundaryPatch::writeDict(Ostream& os) const
{
    this->operator<<(os);
}

Ostream& boundaryPatch::operator<<(Ostream& os) const
{
    os  << name_ << nl << token::BEGIN_BLOCK << nl
        << "    type " << type_ << token::END_STATEMENT << nl
        << "    nFaces " << nFaces_ << token::END_STATEMENT << nl
        << "    startFace " << startFace_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;

    return os;
}

Istream& boundaryPatch::operator>>(Istream& is)
{
    token t;
    is >> name_ >> t;
    is >> t >> type_ >> t;
    is >> t >> nFaces_ >> t;
    is >> t >> startFace_ >> t;
    is >> t;

    return is;
}

void boundaryPatch::operator=(const boundaryPatch& wp)
{
    name_ = wp.name_;
    type_ = wp.type_;
    nFaces_ = wp.nFaces_;
    startFace_ = wp.startFace_;
}

bool boundaryPatch::operator!=(const boundaryPatch& wp) const
{
    if( name_ != wp.name_ )
    {
        return true;
    }
    else if( type_ != wp.name_ )
    {
        return true;
    }
    else if( (nFaces_ != wp.nFaces_) || (startFace_ != wp.startFace_) )
    {
        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
