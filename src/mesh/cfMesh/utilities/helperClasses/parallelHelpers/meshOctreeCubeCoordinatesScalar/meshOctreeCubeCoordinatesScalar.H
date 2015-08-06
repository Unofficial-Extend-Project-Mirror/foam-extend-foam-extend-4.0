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

Class
    meshOctreCubeCoordinatesScalar

Description
    A class containing meshOctreeCubeCoordinates and a scalar value.
    It is used for exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCubeCoordinatesScalar_H
#define meshOctreeCubeCoordinatesScalar_H

#include "scalar.H"
#include "meshOctreeCubeCoordinates.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class meshOctreeCubeCoordinatesScalar Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCubeCoordinatesScalar
{
    // Private data
        //- scalar
        scalar sValue_;

        //- cube coordinates
        meshOctreeCubeCoordinates coordinates_;

    public:

    // Constructors
        //- Null construct
        meshOctreeCubeCoordinatesScalar()
        :
            sValue_(0.0),
            coordinates_()
        {}

        //- Construct from label and cube coordinates
        meshOctreeCubeCoordinatesScalar
        (
            const scalar s,
            const meshOctreeCubeCoordinates& cc
        )
        :
            sValue_(s),
            coordinates_(cc)
        {}

    // Destructor
        ~meshOctreeCubeCoordinatesScalar()
        {}

    // Member functions
        //- return cube label
        inline scalar scalarValue() const
        {
            return sValue_;
        }

        //- return the value
        inline const meshOctreeCubeCoordinates& coordinates() const
        {
            return coordinates_;
        }

    // Member operators

        inline bool operator<(const meshOctreeCubeCoordinatesScalar& scc) const
        {
            if( coordinates_ < scc.coordinates_ )
                return true;

            return false;
        }

        inline void operator=(const meshOctreeCubeCoordinatesScalar& scc)
        {
            sValue_ = scc.sValue_;
            coordinates_ = scc.coordinates_;
        }

        inline bool operator==
        (
            const meshOctreeCubeCoordinatesScalar& scc
        ) const
        {
            if( coordinates_ == scc.coordinates_ )
                return true;

            return false;
        }

        inline bool operator!=
        (
            const meshOctreeCubeCoordinatesScalar& scc
        ) const
        {
            return !this->operator==(scc);
        }

    // Friend operators
        friend Ostream& operator<<
        (
            Ostream& os,
            const meshOctreeCubeCoordinatesScalar& scc
        )
        {
            os << token::BEGIN_LIST;
            os << scc.sValue_ << token::SPACE;
            os << scc.coordinates_ << token::END_LIST;

            // Check state of Ostream
            os.check
            (
                "operator<<(Ostream&, const meshOctreCubeCoordinatesScalar&"
            );

            return os;
        }

        friend Istream& operator>>
        (
            Istream& is,
            meshOctreeCubeCoordinatesScalar& scc
        )
        {
            // Read beginning of meshOctreCubeCoordinatesScalar
            is.readBegin("meshOctreCubeCoordinatesScalar");

            is >> scc.sValue_;
            is >> scc.coordinates_;

            // Read end of meshOctreCubeCoordinatesScalar
            is.readEnd("meshOctreCubeCoordinatesScalar");

            // Check state of Istream
            is.check("operator>>(Istream&, meshOctreCubeCoordinatesScalar");

            return is;
        }
};

//- Specify data associated with meshOctreCubeCoordinatesScalar
//- type is contiguous
template<>
inline bool contiguous<meshOctreeCubeCoordinatesScalar>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
