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
    BlockMatrixCoarsening

Description
    Virtual base class for block matrix AMG coarsening, implemented for
    BlockLduMatrix

Author
    Klas Jareteg, 2012-12-15

SourceFiles
    BlockMatrixCoarsening.C

\*---------------------------------------------------------------------------*/

#ifndef BlockMatrixCoarsening_H
#define BlockMatrixCoarsening_H

#include "runTimeSelectionTables.H"
#include "primitiveFieldsFwd.H"
#include "FieldFields.H"
#include "lduInterfaceFieldPtrsList.H"
#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class coarseBlockAMGLevel;

template<class Type>
class fineBlockAMGLevel;

template<class Type>
class BlockAMGLevel;

/*---------------------------------------------------------------------------*\
                   Class BlockMatrixCoarsening Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockMatrixCoarsening
{
    // Private Data

        //- Dictionary
        const dictionary& dict_;

        //- Group size
        label groupSize_;

        //- Minimum number of coarse level equations
        label minCoarseEqns_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockMatrixCoarsening(const BlockMatrixCoarsening&);

        //- Disallow default bitwise assignment
        void operator=(const BlockMatrixCoarsening&);


public:

    //- Runtime type information
    TypeName("BlockMatrixCoarsening");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockMatrixCoarsening,
            matrix,
            (
                const BlockLduMatrix<Type>& matrix,
                const dictionary& dict,
                const label groupSize,
                const label minCoarseEqns
            ),
            (matrix, dict, groupSize, minCoarseEqns)
        );


    // Selectors

        //- Select given name, group size and matrix
        static autoPtr<BlockMatrixCoarsening<Type> > New
        (
            const word& coarseningType,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const label groupSize,
            const label minCoarseEqns
        );


    // Constructors

        //- Construct from components
        BlockMatrixCoarsening
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const label groupSize,
            const label minCoarseEqns
        )
        :
            dict_(dict),
            groupSize_(groupSize),
            minCoarseEqns_(minCoarseEqns)
        {}


    // Destructor

        virtual ~BlockMatrixCoarsening()
        {}


    // Member Functions

        //- Return access to dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return group size
        label groupSize() const
        {
            return groupSize_;
        }

        //- Return minimum number of coarse level equations
        label minCoarseEqns() const
        {
            return minCoarseEqns_;
        }

        //- Can a coarse level be constructed?
        virtual bool coarsen() const = 0;

        //- Restrict matrix
        virtual autoPtr<BlockAMGLevel<Type> > restrictMatrix
        () const = 0;

        //- Restrict residual
        virtual void restrictResidual
        (
            const Field<Type>& res,
            Field<Type>& coarseRes
        ) const = 0;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            Field<Type>& x,
            const Field<Type>& coarseX
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockMatrixCoarsening.C"
#endif

#define makeBlockMatrixCoarsening(BlockMatrixCoarseningType, typeBlockMatrixCoarseningType) \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typeBlockMatrixCoarseningType, 0);              \
                                                                             \
addToRunTimeSelectionTable(BlockMatrixCoarseningType, typeBlockMatrixCoarseningType, matrix);

#define makeBlockMatrixCoarsenings(blockAmgCoarseningType)                               \
                                                                             \
makeBlockMatrixCoarsening(blockScalarMatrixCoarsening, blockAmgCoarseningType##Scalar);         \
makeBlockMatrixCoarsening(blockVectorMatrixCoarsening, blockAmgCoarseningType##Vector);         \
makeBlockMatrixCoarsening(blockTensorMatrixCoarsening, blockAmgCoarseningType##Tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
