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

#include "meshOctreeCubeBasic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline meshOctreeCubeBasic::meshOctreeCubeBasic()
:
    meshOctreeCubeCoordinates(),
    cubeType_(UNKNOWN),
    procNo_(ALLPROCS)
{}

inline meshOctreeCubeBasic::meshOctreeCubeBasic
(
    const meshOctreeCubeCoordinates& cc
)
:
    meshOctreeCubeCoordinates(cc),
    cubeType_(UNKNOWN),
    procNo_(ALLPROCS)
{}

inline meshOctreeCubeBasic::meshOctreeCubeBasic
(
    const meshOctreeCubeCoordinates& cc,
    const direction cubeType,
    const short procNo
)
:
    meshOctreeCubeCoordinates(cc),
    cubeType_(cubeType),
    procNo_(procNo)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

inline void meshOctreeCubeBasic::setCubeType(const direction ct)
{
    cubeType_ = ct;
}

inline direction meshOctreeCubeBasic::cubeType() const
{
    return cubeType_;
}

inline short meshOctreeCubeBasic::procNo() const
{
    return procNo_;
}

inline void meshOctreeCubeBasic::setProcNo(const short n)
{
    procNo_ = n;
}

inline const meshOctreeCubeCoordinates& meshOctreeCubeBasic::coordinates() const
{
    return *this;
}

inline void meshOctreeCubeBasic::operator=(const meshOctreeCubeBasic& ocb)
{
    cubeType_ = ocb.cubeType_;
    procNo_ = ocb.procNo_;

    static_cast<meshOctreeCubeCoordinates&>(*this) =
        static_cast<const meshOctreeCubeCoordinates&>(ocb);
}

inline Ostream& operator<<
(
    Ostream& os,
    const meshOctreeCubeBasic& ocb
)
{
    os << token::BEGIN_LIST;

    os << ocb.cubeType_ << token::SPACE;
    os << ocb.procNo_ << token::SPACE;
    os << static_cast<const meshOctreeCubeCoordinates&>(ocb);

    os << token::END_LIST;

    // Check state of Ostream
    os.check("operator<<(Ostream&, const meshOctreeCubeBasic");

    return os;
}

inline Istream& operator>>
(
    Istream& is,
    meshOctreeCubeBasic& ocb
)
{
    // Read beginning of meshOctreeCubeBasic
    is.readBegin("meshOctreeCubeBasic");

    label l;
    is >> l;
    ocb.cubeType_ = l;
    is >> l;
    ocb.procNo_ = l;
    is >> static_cast<meshOctreeCubeCoordinates&>(ocb);

    // Read end of meshOctreeCubeBasic
    is.readEnd("meshOctreeCubeBasic");

    // Check state of Istream
    is.check("operator>>(Istream&, meshOctreeCubeBasic");

    return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
