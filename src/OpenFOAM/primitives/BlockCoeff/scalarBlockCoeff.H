/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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

        //- Return as scalar
        const scalar& asScalar() const
        {
            return scalarCoeff_;
        }

        scalar& asScalar()
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
