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
    BlockCoeff specialisation for terminal class.

Description
    Specialisation of a block coefficient for a scalar: always a scalae.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockCoeff_H
#define scalarBlockCoeff_H

#include "BlockCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class BlockCoeff Declaration
\*---------------------------------------------------------------------------*/

template<>
class BlockCoeff<scalar>
:
    public blockCoeffBase
{
public:

    // Public data types

        //- Component type
        typedef scalar Type;

        typedef Type xType;

        //- Coefficient type
        typedef scalar scalarType;


    //- Multiplication trait
    class multiply
    {
    public:

        multiply() {}

        Type operator()(const scalarType& c, const Type& x) const
        {
            return c*x;
        }
        Type operator()(const BlockCoeff<Type>& c, const Type& x) const
        {
            return c.asScalar()*x;
        }

        // Inverse functions
        scalarType inverse(const scalarType& c) const
        {
            return 1.0/c;
        }

        // Triple product of coefficients
        scalarType tripleProduct
        (
            const scalarType& a,
            const scalarType& b,
            const scalarType& c
        ) const
        {
            return a*c/b;
        }
    };


private:

    // Private data

        //- Scalar coefficient
        scalar scalarCoeff_;


public:

    // Constructors

        //- Construct null
        explicit BlockCoeff();

        //- Construct as copy
        BlockCoeff(const BlockCoeff<scalar>&);

        //- Construct from Istream
        BlockCoeff(Istream&);

        //- Clone
        BlockCoeff<scalar> clone() const;


    // Destructor

        ~BlockCoeff();


    // Member functions

        //- Return active type
        blockCoeffBase::activeLevel activeType() const;

        //- Return as scalar
        const scalar& asScalar() const
        {
            return scalarCoeff_;
        }

        scalar& asScalar()
        {
            return scalarCoeff_;
        }

        //- Return as linear
        const scalar& asLinear() const
        {
            return scalarCoeff_;
        }

        scalar& asLinear()
        {
            return scalarCoeff_;
        }

        //- Return as square
        const scalar& asSquare() const
        {
            return scalarCoeff_;
        }

        scalar& asSquare()
        {
            return scalarCoeff_;
        }


        //- Return component
        scalar component(const direction) const;


    // Member operators

        void operator=(const BlockCoeff<scalar>&);


    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const BlockCoeff<scalar>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
