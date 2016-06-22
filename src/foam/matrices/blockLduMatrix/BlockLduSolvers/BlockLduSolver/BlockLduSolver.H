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
    BlockLduSolver

Description
    Linear equation solver for BlockLduMatrix.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockLduSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduSolver_H
#define BlockLduSolver_H

#include "blockLduMatrices.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockSolverPerformance;

/*---------------------------------------------------------------------------*\
                        Class BlockLduSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduSolver
{
    // Private data

        //- Name of field being solved for
        word fieldName_;

        //- Control data dictionary
        dictionary dict_;


protected:

    // Protected data types

    typedef Field<Type> TypeField;


    // Protected data

        //- Matrix
        const BlockLduMatrix<Type>& matrix_;


    // Protected Member Functions

        //- Return dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Read a control parameter from the dictionary
        template<class T>
        void readControl
        (
            const dictionary& dict,
            T& control,
            const word& controlName
        );


public:

    //- Runtime type information
    TypeName("BlockLduSolver");


    // Static data

        //- Large value for the use in solvers
        static const scalar great_;

        //- Small value for the use in solvers
        static const scalar small_;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockLduSolver,
            symMatrix,
            (
                const word& fieldName,
                BlockLduMatrix<Type>& matrix,
                const dictionary& dict
            ),
            (
                fieldName,
                matrix,
                dict
            )
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockLduSolver,
            asymMatrix,
            (
                const word& fieldName,
                BlockLduMatrix<Type>& matrix,
                const dictionary& dict
            ),
            (
                fieldName,
                matrix,
                dict
            )
        );



    // Constructors

        //- Construct from field name and matrix
        BlockLduSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix
        );

        //- Construct from dictionary
        BlockLduSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Selectors

        //- Return a new solver from a dictionary
        static autoPtr<BlockLduSolver<Type> > New
        (
            const word& fieldName,
            BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );

        //- Return a new solver given type
        static autoPtr<BlockLduSolver<Type> > New
        (
            const word& solverName,
            const word& fieldName,
            BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockLduSolver()
        {}


    // Member functions

        //- Return field name
        const word& fieldName() const
        {
            return fieldName_;
        }

        //- Solve
        virtual BlockSolverPerformance<Type> solve
        (
            TypeField& x,
            const TypeField& b
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockLduSolver.C"
#endif


#define makeBlockSolverTypeName(typeSolverType)                              \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typeSolverType, 0);

#define makeBlockSolvers(solverType)                                         \
                                                                             \
makeBlockSolverTypeName(solverType##Scalar);                                 \
makeBlockSolverTypeName(solverType##Vector);                                 \
makeBlockSolverTypeName(solverType##Tensor);


#define addSolverToBlockMatrix(Type, typeSolverType, typeMatrix)             \
                                                                             \
addToRunTimeSelectionTable(block##Type##Solver, typeSolverType, typeMatrix);

#define addSolversToBlockMatrix(solverType, typeMatrix)                      \
                                                                             \
addSolverToBlockMatrix(Scalar, solverType##Scalar, typeMatrix);              \
addSolverToBlockMatrix(Vector, solverType##Vector, typeMatrix);              \
addSolverToBlockMatrix(Tensor, solverType##Tensor, typeMatrix);

#define addSymSolverToBlockMatrix(solverType)                                \
                                                                             \
addSolversToBlockMatrix(solverType, symMatrix);

#define addAsymSolverToBlockMatrix(solverType)                               \
                                                                             \
addSolversToBlockMatrix(solverType, asymMatrix);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
