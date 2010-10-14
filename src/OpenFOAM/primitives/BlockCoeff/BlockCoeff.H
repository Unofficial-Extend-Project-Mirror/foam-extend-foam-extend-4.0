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
    blockCoeff

Description
    Block coefficient combines a scalar, linear and square coefficient
    for different levels of coupling

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef BlockCoeff_H
#define BlockCoeff_H

#include "blockCoeffBase.H"
#include "expandTensor.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockCoeff;

template<class Type>
Ostream& operator<<(Ostream&, const BlockCoeff<Type>&);

/*---------------------------------------------------------------------------*\
                          Class BlockCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockCoeff
:
    public blockCoeffBase
{
public:

    // Public data types

        //- Component type
        typedef Type xType;
        typedef Field<xType> xTypeField;

        //- Coefficient type
        typedef typename pTraits<Type>::cmptType scalarType;
        typedef Type linearType;
        typedef typename outerProduct<Type, Type>::type squareType;

        //- Field type
        typedef Field<scalarType> scalarTypeField;
        typedef Field<linearType> linearTypeField;
        typedef Field<squareType> squareTypeField;


    //- Multiplication trait
    class multiply
    {
    public:

        multiply() {}

        // Coefficient times type multiplication

            Type operator()(const scalarType& c, const Type& x) const
            {
                return c*x;
            }

            Type operator()(const linearType& c, const Type& x) const
            {
                return cmptMultiply(c, x);
            }

            Type operator()(const squareType& c, const Type& x) const
            {
                return (c & x);
            }

            Type operator()(const BlockCoeff<Type>& c, const Type& x) const
            {
                if (c.scalarCoeffPtr_)
                {
                    return operator()(*c.scalarCoeffPtr_, x);
                }
                else if (c.linearCoeffPtr_)
                {
                    return operator()(*c.linearCoeffPtr_, x);
                }
                else if (c.squareCoeffPtr_)
                {
                    return operator()(*c.squareCoeffPtr_, x);
                }
                else
                {
                    return pTraits<Type>::zero;
                }
            }


        // Transpose functions

            scalarType transpose(const scalarType& c) const
            {
                return c;
            }

            linearType transpose(const linearType& c) const
            {
                return c;
            }

            squareType transpose(const squareType& c) const
            {
                return c.T();
            }


        // Inverse functions

            scalarType inverse(const scalarType& c) const
            {
                return 1.0/c;
            }

            linearType inverse(const linearType& c) const
            {
                return cmptDivide(pTraits<linearType>::one, c);
            }

            squareType inverse(const squareType& c) const
            {
                return inv(c);
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

            linearType tripleProduct
            (
                const scalarType& a,
                const linearType& b,
                const scalarType& c
            ) const
            {
                return a*c*inverse(b);
            }

            linearType tripleProduct
            (
                const linearType& a,
                const linearType& b,
                const linearType& c
            ) const
            {
                return cmptDivide(cmptMultiply(a, c), b);
            }

            squareType tripleProduct
            (
                const scalarType& a,
                const squareType& b,
                const scalarType& c
            ) const
            {
                return a*c*inv(b);
            }

            squareType tripleProduct
            (
                const linearType& a,
                const squareType& b,
                const linearType& c
            ) const
            {
                squareType result;
                linearType sac = cmptMultiply(a, c);

                expandLinear(result, sac);
                return result & inv(b);
            }

            squareType tripleProduct
            (
                const squareType& a,
                const squareType& b,
                const squareType& c
            ) const
            {
                return (a & inv(b)) & c;
            }
    };


private:

    // Private data

        //- Scalar coefficient
        mutable scalarType* scalarCoeffPtr_;

        //- Linear coefficient
        mutable linearType* linearCoeffPtr_;

        //- Square coefficient
        mutable squareType* squareCoeffPtr_;


    // Private Member Functions

        //- Promote to scalar
        scalarType& toScalar();

        //- Promote to linear
        linearType& toLinear();

        //- Promote to square
        squareType& toSquare();


public:

    // Constructors

        //- Construct null
        explicit BlockCoeff();

        //- Construct as copy
        BlockCoeff(const BlockCoeff<Type>&);

        //- Construct from Istream
        BlockCoeff(Istream&);

        //- Clone
        BlockCoeff<Type> clone() const;


    // Destructor

        ~BlockCoeff();

        //- Clear data
        void clear();


    // Member functions

        //- Return active type
        blockCoeffBase::activeLevel activeType() const;

        //- Check pointers: only one type should be active (debug only)
        void checkActive() const;

        // Return as typed.  Fails when asked for the incorrect type

            //- Return as scalar
            const scalarType& asScalar() const;
            scalarType& asScalar();

            //- Return as linear
            const linearType& asLinear() const;
            linearType& asLinear();

            //- Return as square
            const squareType& asSquare() const;
            squareType& asSquare();


        //- Return component
        scalarType component(const direction) const;


    // Member operators

        void operator=(const BlockCoeff<Type>&);


    // IOstream operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const BlockCoeff<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
