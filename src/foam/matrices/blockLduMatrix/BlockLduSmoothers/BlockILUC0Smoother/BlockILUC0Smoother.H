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
    BlockILUC0Smoother

Description
    ILU Cp smoother

SourceFiles
    blockILUC0Smoothers.C

\*---------------------------------------------------------------------------*/

#ifndef BlockILUC0Smoother_H
#define BlockILUC0Smoother_H

#include "BlockLduSmoother.H"
#include "blockILUC0Precons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BlockILUC0Smoother Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockILUC0Smoother
:
    public BlockLduSmoother<Type>
{
    // Private Data

        //- ILUCp preconditioner
        BlockILUC0Precon<Type> precon_;

        //- Correction array
        mutable Field<Type> xCorr_;

        //- Residual array
        mutable Field<Type> residual_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockILUC0Smoother(const BlockILUC0Smoother&);

        //- Disallow default bitwise assignment
        void operator=(const BlockILUC0Smoother&);


public:

    //- Runtime type information
    TypeName("ILUC0");


    // Constructors

        //- Construct from components
        BlockILUC0Smoother
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        )
        :
            BlockLduSmoother<Type>(matrix),
            precon_(matrix, dict),
            xCorr_(matrix.lduAddr().size()),
            residual_(matrix.lduAddr().size())
        {}


    //- Destructor
    virtual ~BlockILUC0Smoother()
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
