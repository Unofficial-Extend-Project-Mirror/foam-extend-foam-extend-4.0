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
    BlockCoeff

Description
    Template for the terminal decoupled class.  It is designed to avoid
    endless expansion of tensor order by excluding block coupling at the
    terminal type level.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef DecoupledBlockCoeff_H
#define DecoupledBlockCoeff_H

#include "blockCoeffBase.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class DecoupledBlockCoeff;

template<class Type>
Ostream& operator<<(Ostream&, const DecoupledBlockCoeff<Type>&);

/*---------------------------------------------------------------------------*\
                      Class DecoupledBlockCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DecoupledBlockCoeff
:
    public blockCoeffBase
{
public:

    // Public data types

        //- Component type
        typedef Type xType;
        typedef Field<xType> xTypeField;

        //- Coefficient type
        typedef scalar scalarType;
        typedef Type linearType;

        //- Field type
        typedef Field<scalarType> scalarTypeField;
        typedef Field<linearType> linearTypeField;


    //- Multiplication trait
    class multiply
    {
    public:

        multiply() {}

        Type operator()(const scalarType& c, const Type& x) const
        {
            return c*x;
        }

        Type operator()(const linearType& c, const Type& x) const
        {
            return cmptMultiply(c, x);
        }

        Type operator()(const DecoupledBlockCoeff<Type>& c, const Type& x) const
        {
            if (c.scalarCoeffPtr_)
            {
                return operator()(*c.scalarCoeffPtr_, x);
            }
            else if (c.linearCoeffPtr_)
            {
                return operator()(*c.linearCoeffPtr_, x);
            }
            else
            {
                return pTraits<Type>::zero;
            }
        }


        // Coefficient times coefficient multiplication. Needed for BlockILUCp
        // preconditioner. VV, 12/Jul/2015.

            scalarType activeTypeMultiply
            (
                const scalarType& a,
                const scalarType& b
            ) const
            {
                return a*b;
            }

            linearType activeTypeMultiply
            (
                const linearType& a,
                const linearType& b
            ) const
            {
                return cmptMultiply(a, b);
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
                const linearType& a,
                const linearType& b,
                const linearType& c
            ) const
            {
                return cmptDivide(cmptMultiply(a, c), b);
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
    };


private:

    // Private data

        //- Scalar coefficient
        mutable scalarType* scalarCoeffPtr_;

        //- Linear coefficient
        mutable linearType* linearCoeffPtr_;


    // Private Member Functions

        //- Promote to scalar
        scalarType& toScalar();

        //- Promote to linear
        linearType& toLinear();


public:

    // Constructors

        //- Construct null
        explicit DecoupledBlockCoeff();

        //- Construct as copy
        DecoupledBlockCoeff(const DecoupledBlockCoeff<Type>&);

        //- Construct from Istream
        DecoupledBlockCoeff(Istream&);

        //- Clone
        DecoupledBlockCoeff<Type> clone() const;


    // Destructor

        ~DecoupledBlockCoeff();

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


        //- Return component
        scalarType component(const direction) const;


    // Member operators

        void operator=(const DecoupledBlockCoeff<Type>&);


    // IOstream operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const DecoupledBlockCoeff<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DecoupledBlockCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
