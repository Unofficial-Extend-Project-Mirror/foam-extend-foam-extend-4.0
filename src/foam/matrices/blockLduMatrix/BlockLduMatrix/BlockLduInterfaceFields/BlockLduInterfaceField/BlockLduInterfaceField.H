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
    BlockLduInterfaceField

Description
    Abstract base class for interface fields with block coefficients

Author
    Ivor Clifford, ivor.clifford@gmail.com

SourceFiles
    BlockLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterfaceField_H
#define BlockLduInterfaceField_H

#include "lduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

template<class Type>
class CoeffField;

/*---------------------------------------------------------------------------*\
                    Class BlockLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterfaceField
:
    public lduInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduInterfaceField(const BlockLduInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterfaceField&);


public:

    //- Runtime Type information
    TypeName("BlockLduInterfaceField");

    // Constructors

        //- Construct given coupled patch
        BlockLduInterfaceField(const lduInterface& patch)
        :
            lduInterfaceField(patch)
        {}


    //- Destructor
    virtual ~BlockLduInterfaceField()
    {}


    // Member Functions

        // Coupled interface matrix update, scalar matrices
        // Used to override base class lduInterfaceField virtual functions

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const = 0;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const
            {}


        // Block coupled interface matrix update

            //- Transform neighbour field
            virtual void transformCoupleField
            (
                Field<Type>& pnf
            ) const = 0;

            // KRJ: 2013-02-08: Must be virtual!
            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const = 0;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
