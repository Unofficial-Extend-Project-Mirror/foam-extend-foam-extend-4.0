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

Class
    Foam::Vector2D

Description
    Templated 2D Vector derived from VectorSpace adding construction from
    2 components, element access using x() and y() member functions and
    the inner-product (dot-product).

SourceFiles
    Vector2DTemplateI.H

\*---------------------------------------------------------------------------*/

#ifndef Vector2DTemplate_H
#define Vector2DTemplate_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Vector2D Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt>
class Vector2D
:
    public VectorSpace<Vector2D<Cmpt>, Cmpt, 2>
{

public:

    // Member constants

        enum
        {
            rank = 1 // Rank of Vector2D is 1
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const Vector2D zero;
        static const Vector2D one;
        static const Vector2D max;
        static const Vector2D min;


    //- Component labeling enumeration
    enum components { X, Y };


    // Constructors

        //- Construct null
        inline Vector2D();

        //- Construct given VectorSpace
        inline Vector2D(const VectorSpace<Vector2D<Cmpt>, Cmpt, 2>&);

        //- Construct given two components
        inline Vector2D(const Cmpt& vx, const Cmpt& vy);

        //- Construct from Istream
        inline Vector2D(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& x() const;
            inline const Cmpt& y() const;

            inline Cmpt& x();
            inline Cmpt& y();

            //- Return i-th component.  Consistency with VectorN
            inline const Cmpt& operator()
            (
                const direction i
            ) const;

            //- Return i-th component.  Consistency with VectorN
            inline  Cmpt& operator()
            (
                const direction i
            );


        // Operators

            //- perp dot product (dot product with perpendicular vector)
            inline scalar perp(const Vector2D<Cmpt>& b) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "Vector2DTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
