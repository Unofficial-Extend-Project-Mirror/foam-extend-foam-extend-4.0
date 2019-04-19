/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    BlockILUSmoother

Description
    ILU smoother

SourceFiles
    blockILUSmoothers.C

\*---------------------------------------------------------------------------*/

#ifndef BlockILUSmoother_H
#define BlockILUSmoother_H

#include "BlockLduSmoother.H"
#include "blockCholeskyPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BlockILUSmoother Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockILUSmoother
:
    public BlockLduSmoother<Type>
{
    // Private Data

        //- Cholesky preconditioner
        BlockCholeskyPrecon<Type> precon_;

        //- Correction array
        mutable Field<Type> xCorr_;

        //- Residual array
        mutable Field<Type> residual_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockILUSmoother(const BlockILUSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const BlockILUSmoother&);


public:

    //- Runtime type information
    TypeName("ILU");


    // Constructors

        //- Construct from components
        BlockILUSmoother
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        )
        :
            BlockLduSmoother<Type>(matrix),
            precon_(matrix),
            xCorr_(matrix.lduAddr().size()),
            residual_(matrix.lduAddr().size())
        {}


    //- Destructor
    virtual ~BlockILUSmoother()
    {}


    // Member Functions

        //- Execute smoothing
        virtual void smooth
        (
            Field<Type>& x,
            const Field<Type>& b,
            const label nSweeps
        ) const
        {
            for (label sweep = 0; sweep < nSweeps; sweep++)
            {
                // Calculate residual
                this-> matrix_.Amul(residual_, x);

                // residual = b - Ax
                forAll (b, i)
                {
                    residual_[i] = b[i] - residual_[i];
                }

                // Note.
                // For the ILU smoother it is not necessary to initialise
                // xCorr to zero, because it gets over-written by the product
                // but for more complex smoothers it is
                // HJ, 15/Jun/2017
                // xCorr_ = pTraits<Type>::zero;
                precon_.precondition(xCorr_, residual_);

                // Add correction to x
                x += xCorr_;
            }
        }

        //- Re-initialise preconditioner after matrix coefficient update
        virtual void initMatrix()
        {
            precon_.initMatrix();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
