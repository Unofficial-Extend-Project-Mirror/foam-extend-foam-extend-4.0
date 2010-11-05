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
    CoeffField<T>

Description
    Generic coefficient field type.  Used in BlockLduMatrix.  HJ, 2/Apr/2005

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    CoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef CoeffField_H
#define CoeffField_H

#include "VectorSpace.H"
#include "primitiveFields.H"
#include "blockCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class CoeffField;

template<class Type>
Ostream& operator<<(Ostream&, const CoeffField<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const tmp<CoeffField<Type> >&);


/*---------------------------------------------------------------------------*\
                          Class CoeffField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CoeffField
:
    public refCount
{
public:

    // Public data types

        //- Component type
        typedef typename BlockCoeff<Type>::scalarType scalarType;
        typedef typename BlockCoeff<Type>::linearType linearType;
        typedef typename BlockCoeff<Type>::squareType squareType;

        //- Field type
        typedef typename BlockCoeff<Type>::scalarTypeField scalarTypeField;
        typedef typename BlockCoeff<Type>::linearTypeField linearTypeField;
        typedef typename BlockCoeff<Type>::squareTypeField squareTypeField;


private:

    // Private data

        //- Scalar coefficient
        mutable scalarTypeField* scalarCoeffPtr_;

        //- Linear coefficient
        mutable linearTypeField* linearCoeffPtr_;

        //- Square coefficient
        mutable squareTypeField* squareCoeffPtr_;

        //- Number of elements
        label size_;


    // Private Member Functions

        //- Check size for arithmetic operations: resizing is not allowed!
        template<class Type2>
        inline void checkSize(const UList<Type2>&) const;

        //- Promote to scalar
        scalarTypeField& toScalar();

        //- Promote to linear
        linearTypeField& toLinear();

        //- Promote to square
        squareTypeField& toSquare();


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct given size
        explicit CoeffField(const label);

        //- Construct as copy
        CoeffField(const CoeffField<Type>&);

        //- Construct from Istream
        explicit CoeffField(Istream&);

        //- Clone
        tmp<CoeffField<Type> > clone() const;


    // Destructor

        ~CoeffField();

        //- Clear data
        void clear();


    // Member functions

        //- Return size
        inline label size() const;

        //- Return active type
        blockCoeffBase::activeLevel activeType() const;

        //- Check pointers: only one type should be active (debug only)
        void checkActive() const;

        //- Negate this field
        void negate();

        //- Return the field transpose
        tmp<CoeffField<Type> > transpose() const;


        // Return as typed.  Fails when asked for the incorrect type

            //- Return as scalar field
            const scalarTypeField& asScalar() const;

            //- Return as linear field
            const linearTypeField& asLinear() const;

            //- Return as square field
            const squareTypeField& asSquare() const;


        // Return as typed.  Fails when asked for demotion

            //- Return as scalar field
            scalarTypeField& asScalar();

            //- Return as linear field
            linearTypeField& asLinear();

            //- Return as square field
            squareTypeField& asSquare();


        //- Return component
        tmp<scalarTypeField> component(const direction) const;

        //- Return coefficient as block
        BlockCoeff<Type> getCoeff(const label index) const;

        //- Set coefficient from a block
        void setCoeff(const label index, const BlockCoeff<Type>& coeff);


        // Subset operations

            //- Get subset with offset and size and store in given field
            void getSubset
            (
                CoeffField<Type>& f,
                const label start,
                const label size
            ) const;

            //- Get subset with addressing and store in given field
            void getSubset
            (
                CoeffField<Type>& f,
                const labelList& addr
            ) const;


            //- Set subset with offset and size from given field
            void setSubset
            (
                const CoeffField<Type>& f,
                const label start,
                const label size
            );

            //- Get subset with addressing and store in target field
            void setSubset
            (
                const CoeffField<Type>& f,
                const labelList& addr
            );

            //- Zero out subset with offset and size
            void zeroOutSubset
            (
                const label start,
                const label size
            );

            //- Zero out subset with addressing
            void zeroOutSubset
            (
                const labelList& addr
            );

            //- Add subset with addressing to field
            void addSubset
            (
                const CoeffField<Type>& f,
                const labelList& addr
            );

            //- Subtract subset with addressing to field
            void subtractSubset
            (
                const CoeffField<Type>& f,
                const labelList& addr
            );
    // Member operators

        void operator=(const CoeffField<Type>&);
        void operator=(const tmp<CoeffField<Type> >&);

        void operator=(const scalarTypeField&);
        void operator=(const tmp<scalarTypeField>&);
        void operator=(const linearTypeField&);
        void operator=(const tmp<linearTypeField>&);
        void operator=(const squareTypeField&);
        void operator=(const tmp<squareTypeField>&);


        void operator+=(const CoeffField<Type>&);
        void operator+=(const tmp<CoeffField<Type> >&);

        void operator+=(const scalarTypeField&);
        void operator+=(const tmp<scalarTypeField>&);
        void operator+=(const linearTypeField&);
        void operator+=(const tmp<linearTypeField>&);
        void operator+=(const squareTypeField&);
        void operator+=(const tmp<squareTypeField>&);

        void operator-=(const CoeffField<Type>&);
        void operator-=(const tmp<CoeffField<Type> >&);

        void operator-=(const scalarTypeField&);
        void operator-=(const tmp<scalarTypeField>&);
        void operator-=(const linearTypeField&);
        void operator-=(const tmp<linearTypeField>&);
        void operator-=(const squareTypeField&);
        void operator-=(const tmp<squareTypeField>&);

        void operator*=(const UList<scalar>&);
        void operator*=(const tmp<Field<scalar> >&);
        void operator*=(const scalar&);

        void operator/=(const UList<scalar>&);
        void operator/=(const tmp<Field<scalar> >&);
        void operator/=(const scalar&);

        tmp<CoeffField<Type> > operator-();

    // IOstream operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const CoeffField<Type>&
        );

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const tmp<CoeffField<Type> >&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoeffFieldFunctions.H"

#ifdef NoRepository
#   include "CoeffField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
