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

#include "DynList.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Construct from Istream
template<class T, Foam::label staticSize>
Foam::DynList<T, staticSize>::DynList(Istream&)
:
    dataPtr_(nullptr),
    nAllocated_(0),
    staticData_(),
    nextFree_(0)
{
    FatalErrorIn
    (
        "template<class T, Foam::label staticSize>"
        "\nFoam::DynList<T, staticSize>::DynList(Istream& is)"
    ) << "Not implemented" << exit(FatalError);
}


template<class T, Foam::label staticSize>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::DynList<T, staticSize>& DL
)
{
    UList<T> helper(DL.dataPtr_, DL.nextFree_);
    os << helper;

    return os;
}


template<class T, Foam::label staticSize>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::DynList<T, staticSize>& DL
)
{
    FatalErrorIn
    (
        "template<class T, Foam::label staticSize>"
        "\nFoam::Istream& Foam::operator>>"
        "(Foam::Istream& is, Foam::DynList<T, staticSize>& DL)"
    ) << "Not implemented" << exit(FatalError);

    UList<T> helper(DL.dataPtr_, DL.nextFree_);
    //is >> static_cast<List<T>&>(DL);
    is >> helper;
    DL.nextFree_ = helper.size();

    return is;
}


// ************************************************************************* //
