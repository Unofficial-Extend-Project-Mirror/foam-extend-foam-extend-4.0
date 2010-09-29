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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    BlockConstraint

Description
    A storage mechanism which allows setting of the fixed value and
    consequently recovering the equation for a single row of the matrix as
    well as the b. The equation is taken out of the matrix using a
    variant of compact matrix storage mechanism.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef BlockConstraint_H
#define BlockConstraint_H

#include "BlockCoeff.H"
#include "coeffFields.H"
#include "blockLduMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class BlockConstraint;

template<class Type>
Ostream& operator<<(Ostream&, const BlockConstraint<Type>&);


/*---------------------------------------------------------------------------*\
                        Class BlockConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockConstraint
{
    // Public data types

    typedef CoeffField<Type> TypeCoeffField;
    typedef Field<Type> TypeField;

    // Private data

        //- Matrix row ID
        label rowID_;

        //- Fixed value
        Type value_;

        //- Fixed components (0-1) 1 = fixed, 0 = free
        Type fixedComponents_;

        //- Are matrix coefficients set?
        bool matrixCoeffsSet_;

        //- Diagonal coefficient
        BlockCoeff<Type> diagCoeff_;

        //- Right-hand side
        Type b_;

        //- Upper coefficients, owner side
        TypeCoeffField* upperCoeffsOwnerPtr_;

        //- Upper coefficients, neighbour side
        TypeCoeffField* upperCoeffsNeighbourPtr_;

        //- Lower coefficients, owner side
        TypeCoeffField* lowerCoeffsOwnerPtr_;

        //- Lower coefficients, neighbour side
        TypeCoeffField* lowerCoeffsNeighbourPtr_;


public:


    // Constructors

        //- Construct from components
        BlockConstraint
        (
            const label row,
            const Type value,
            const Type& fixedCmpts = pTraits<Type>::one
        );

        //- Construct as copy
        BlockConstraint(const BlockConstraint&);

        //- Construct from Istream
        BlockConstraint(Istream&);


    // Destructor

        ~BlockConstraint();


    // Member Functions

        //- Return matrix row ID
        label rowID() const
        {
            return rowID_;
        }

        //- Return fixed value
        Type value() const
        {
            return value_;
        }

        //- Return map of fixed components
        const Type& fixedComponents() const
        {
            return fixedComponents_;
        }

        //- Return diagonal coefficient
        const BlockCoeff<Type>& diagCoeff() const;

        //- Return right-hand side
        const Type& b() const;

        //- Return off-diagonal coefficients
        const TypeCoeffField& upperCoeffsOwner() const;

        const TypeCoeffField& upperCoeffsNeighbour() const;

        const TypeCoeffField& lowerCoeffsOwner() const;

        const TypeCoeffField& lowerCoeffsNeighbour() const;

        //- Combine with existing equation
        void combine(const BlockConstraint<Type>&);

        //- Set matrix coefficients
        void setMatrix
        (
            const BlockLduMatrix<Type>& matrix,
            const TypeField& x,
            const TypeField& b
        );


        //- Eliminate equation
        void eliminateEquation
        (
            BlockLduMatrix<Type>& matrix,
            TypeField& b
        ) const;

        //- Set x, b and diagonal in eliminated equation
        void setSourceDiag
        (
            BlockLduMatrix<Type>&,
            Field<Type>& x,
            Field<Type>& b
        ) const;


        //- Reconstruct matrix coefficients
        void reconstructMatrix(BlockLduMatrix<Type>&) const;

        //- Clear matrix coefficients
        void clearMatrix();


    // Member Operators

        void operator=(const BlockConstraint<Type>&);


    // Ostream Operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const BlockConstraint<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "BlockConstraint.C"
#    include "BlockConstraintTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
