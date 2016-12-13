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
    DiagTensorN

Description
    Templated NXN DiagTensor derived from VectorSpace adding construction from
    N components, and the inner-product (dot-product) and outer-product operators.

SourceFiles
    DiagTensorNI.H

\*---------------------------------------------------------------------------*/

#ifndef DiagTensorN_H
#define DiagTensorN_H

#include "VectorN.H"
#include "SphericalTensorN.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DiagTensor Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt, int length>
class DiagTensorN
:
    public VectorSpace<DiagTensorN<Cmpt, length>, Cmpt, length>
{

public:

    // Member constants

        enum
        {
            rank = 2, // Rank of DiagTensor is 2
            rowLength = length
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const DiagTensorN zero;
        static const DiagTensorN one;
        static const DiagTensorN I;


    // Constructors

        //- Construct null
        inline DiagTensorN();

        //- Construct given VectorSpace
        inline DiagTensorN(const VectorSpace<DiagTensorN<Cmpt, length>, Cmpt, length>&);

        //- Construct from Istream
        inline DiagTensorN(Istream&);

        //- Construct given component value.  Special use only!
        explicit inline DiagTensorN(const Cmpt& tx);


    // Member Functions

        //- Diagonal
        inline DiagTensorN<Cmpt, length> diag() const;

        //- Transpose
        inline DiagTensorN<Cmpt, length> T() const;

    // Member Operators

        //- Assign to a SphericalTensorN
        inline void operator=(const SphericalTensorN<Cmpt, length>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "DiagTensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
