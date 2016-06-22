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
    VectorN

Description
    Templated 2D Vector derived from VectorSpace adding construction from
    2 components, element access using x() and y() member functions and
    the inner-product (dot-product).

SourceFiles
    VectorNI.H

\*---------------------------------------------------------------------------*/

#ifndef VectorN_H
#define VectorN_H

#include "VectorSpace.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class VectorN Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt, int length>
class VectorN
:
    public VectorSpace<VectorN<Cmpt, length>, Cmpt, length>
{

public:

    // Member constants

        enum
        {
            rank = 1 // Rank of VectorN is 1
        };


    // Static data members

        static const char* const typeName;
        static const VectorN zero;
        static const VectorN one;
        static const VectorN max;
        static const VectorN min;


    // Constructors

        //- Construct null
        inline VectorN();

        //- Construct given VectorSpace
        inline VectorN
        (
            const VectorSpace<VectorN<Cmpt, length>, Cmpt, length>&
        );

        //- Construct given component value.  Special use only!
        explicit inline VectorN(const Cmpt& vx);

        //- Construct from Istream
        inline VectorN(Istream&);


    // Member Functions

        //- Return ith component
        inline const Cmpt& operator()(const direction i) const;

        //- Return access to ith component
        inline Cmpt& operator()(const direction i);

        //- Component-wise multiply
        inline VectorN<Cmpt, length> cmptMultiply
        (
            const VectorN<Cmpt, length>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "VectorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
