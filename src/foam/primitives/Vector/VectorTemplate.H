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
    Foam::Vector

Description
    Templated 3D Vector derived from VectorSpace adding construction from
    3 components, element access using x(), y() and z() member functions and
    the inner-product (dot-product) and cross product operators.

    A centre() member function which returns the Vector for which it is called
    is defined so that point which is a typedef to Vector\<scalar\> behaves as
    other shapes in the shape hierachy.

SourceFiles
    VectorTemplateI.H

\*---------------------------------------------------------------------------*/

#ifndef VectorTemplate_H
#define VectorTemplate_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T> class List;

/*---------------------------------------------------------------------------*\
                           Class Vector Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt>
class Vector
:
    public VectorSpace<Vector<Cmpt>, Cmpt, 3>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Vector<label> labelType;


    // Member constants

        enum
        {
            rank = 1 // Rank of Vector is 1
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const Vector zero;
        static const Vector one;
        static const Vector max;
        static const Vector min;


    //- Component labeling enumeration
    enum components { X, Y, Z };


    // Constructors

        //- Construct null
        inline Vector();

        //- Construct given VectorSpace
        inline Vector(const VectorSpace<Vector<Cmpt>, Cmpt, 3>&);

        //- Construct given three components
        inline Vector(const Cmpt& vx, const Cmpt& vy, const Cmpt& vz);

        //- Construct from Istream
        inline Vector(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& x() const;
            inline const Cmpt& y() const;
            inline const Cmpt& z() const;

            inline Cmpt& x();
            inline Cmpt& y();
            inline Cmpt& z();

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


        //- Return *this (used for point which is a typedef to Vector<scalar>.
        inline const Vector<Cmpt>& centre
        (
            const Foam::List<Vector<Cmpt> >&
        ) const;
};


template<class Cmpt>
class typeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VectorTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
