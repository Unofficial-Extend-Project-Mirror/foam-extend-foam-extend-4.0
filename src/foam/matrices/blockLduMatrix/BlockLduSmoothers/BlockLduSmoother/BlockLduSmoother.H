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
    BlockLduSmoother

Description
    Block LDU matrix smoother virtual base class

SourceFiles
    newBlockLduSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduSmoother_H
#define BlockLduSmoother_H

#include "blockLduMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BlockLduSmoother Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduSmoother
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduSmoother(const BlockLduSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduSmoother&);


protected:

    // Protected data

        //- Matrix reference
        const BlockLduMatrix<Type>& matrix_;


public:

    //- Runtime type information
    TypeName("BlockLduSmoother");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockLduSmoother,
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
        BlockLduSmoother(const BlockLduMatrix<Type>& matrix)
        :
            matrix_(matrix)
        {}


    // Selectors

        //- Select given matrix and dictionary
        static autoPtr<BlockLduSmoother<Type> > New
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const word keyword = word("smoother")
        );


    // Destructor

        virtual ~BlockLduSmoother()
        {}


    // Member Functions

        //- Execute smoothing
        virtual void smooth
        (
            Field<Type>& x,
            const Field<Type>& b,
            const label nSweeps
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockLduSmoother.C"
#endif

#define makeBlockSmoother(SmootherType, typeSmootherType)                    \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typeSmootherType, 0);                    \
                                                                             \
addToRunTimeSelectionTable(SmootherType, typeSmootherType, dictionary);

#define makeBlockSmoothers(smootherType)                                       \
                                                                             \
makeBlockSmoother(blockScalarSmoother, smootherType##Scalar);                  \
makeBlockSmoother(blockVectorSmoother, smootherType##Vector);                  \
makeBlockSmoother(blockTensorSmoother, smootherType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
