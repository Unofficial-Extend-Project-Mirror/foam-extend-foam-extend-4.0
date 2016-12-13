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
    BlockILUCpSmoother

Description
    ILU Cp smoother

SourceFiles
    blockILUCpSmoothers.C

\*---------------------------------------------------------------------------*/

#ifndef BlockILUCpSmoother_H
#define BlockILUCpSmoother_H

#include "BlockLduSmoother.H"
#include "blockILUCpPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BlockILUCpSmoother Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockILUCpSmoother
:
    public BlockLduSmoother<Type>
{
    // Private Data

        //- ILUCp preconditioner
        BlockILUCpPrecon<Type> precon_;

        //- Correction array
        mutable Field<Type> xCorr_;

        //- Residual array
        mutable Field<Type> residual_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockILUCpSmoother(const BlockILUCpSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const BlockILUCpSmoother&);


public:

    //- Runtime type information
    TypeName("ILUCp");


    // Constructors

        //- Construct from components
        BlockILUCpSmoother
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
    virtual ~BlockILUCpSmoother()
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
