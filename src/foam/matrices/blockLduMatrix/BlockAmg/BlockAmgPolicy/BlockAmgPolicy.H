/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    BlockAmgPolicy

Description
    Virtual base class for AMG coarsening policy, implemented for
    BlockLduMatrix

Author
    Klas Jareteg, 2012-12-15

SourceFiles
    BlockAmgPolicy.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAmgPolicy_H
#define BlockAmgPolicy_H

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
class coarseBlockAmgLevel;

template<class Type>
class fineBlockAmgLevel;

template<class Type>
class BlockAmgLevel;

/*---------------------------------------------------------------------------*\
                           Class BlockAmgPolicy Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAmgPolicy
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
        BlockAmgPolicy(const BlockAmgPolicy&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAmgPolicy&);


public:

        //- Runtime type information
        TypeName("BlockAmgPolicy");

        // KRJ: 2012-12-15: Replaced by the above function
        //- Runtime type information
        //virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                BlockAmgPolicy,
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
        static autoPtr<BlockAmgPolicy<Type> > New
        (
            const word& policyType,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const label groupSize,
            const label minCoarseEqns
        );


    // Constructors

        //- Construct from components
        BlockAmgPolicy
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

        virtual ~BlockAmgPolicy()
        {}


    // Member Functions

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
        virtual autoPtr<BlockLduMatrix<Type> > restrictMatrix
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
#   include "BlockAmgPolicy.C"
#endif

#define makeBlockAmgPolicy(BlockAmgPolicyType, typeBlockAmgPolicyType) \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typeBlockAmgPolicyType, 0);              \
                                                                             \
addToRunTimeSelectionTable(BlockAmgPolicyType, typeBlockAmgPolicyType, matrix);

#define makeBlockAmgPolicys(blockAmgPolicyType)                               \
                                                                             \
makeBlockAmgPolicy(blockScalarAmgPolicy, blockAmgPolicyType##Scalar);         \
makeBlockAmgPolicy(blockVectorAmgPolicy, blockAmgPolicyType##Vector);         \
makeBlockAmgPolicy(blockTensorAmgPolicy, blockAmgPolicyType##Tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
