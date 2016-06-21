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
    BlockLduPrecon

Description
    Block LDU matrix preconditioner virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    newBlockLduPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduPrecon_H
#define BlockLduPrecon_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockLduPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduPrecon
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduPrecon(const BlockLduPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduPrecon&);


protected:

    // Protected data

        //- Matrix reference
        const BlockLduMatrix<Type>& matrix_;


public:

    //- Runtime type information
    TypeName("BlockLduPrecon");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockLduPrecon,
            dictionary,
            (
                const BlockLduMatrix<Type>& matrix,
                const dictionary& dict
            ),
            (
                matrix,
                dict
            )
        );


    // Constructors

        //- Construct from matrix
        BlockLduPrecon(const BlockLduMatrix<Type>& matrix)
        :
            matrix_(matrix)
        {}


    // Selectors

        //- Select given matrix and dictionary
        static autoPtr<BlockLduPrecon<Type> > New
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const word keyword = word("preconditioner")
        );


    // Destructor

        virtual ~BlockLduPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const = 0;

        //- Execute preconditioning on a transposed matrix
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const
        {
            notImplemented
            (
                type() +"::preconditionT"
                "(Field<Type>& xT, const Field<Type>& bT) const"
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockLduPrecon.C"
#endif

#define makeBlockPrecon(PreconType, typePreconType)                          \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typePreconType, 0);                      \
                                                                             \
addToRunTimeSelectionTable(PreconType, typePreconType, dictionary);

#define makeBlockPrecons(preconType)                                         \
                                                                             \
makeBlockPrecon(blockScalarPrecon, preconType##Scalar);                      \
makeBlockPrecon(blockVectorPrecon, preconType##Vector);                      \
makeBlockPrecon(blockTensorPrecon, preconType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
