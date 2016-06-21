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
    fineBlockAMGLevel

Description
    Fine AMG level refers to matrix, x and b data from the solver or
    preconditioner, implemented for the BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

SourceFiles
    fineBlockAMGLevel.C

\*---------------------------------------------------------------------------*/

#ifndef fineBlockAMGLevel_H
#define fineBlockAMGLevel_H

#include "BlockAMGLevel.H"
#include "BlockLduSmoother.H"
#include "BlockLduMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class coarseBlockAMGLevel;

template<class Type>
class BlockMatrixCoarsening;

/*---------------------------------------------------------------------------*\
                        Class fineBlockAMGLevel Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fineBlockAMGLevel
:
    public BlockAMGLevel<Type>
{
    // Private data

         //- Reference to matrix
        const BlockLduMatrix<Type>& matrix_;

        //- Dictionary
        const dictionary& dict_;

        //- AMG coarsening
        autoPtr<BlockMatrixCoarsening<Type> > coarseningPtr_;

        //- Smoother
        autoPtr<BlockLduSmoother<Type> > smootherPtr_;

        //- Ax buffer
        mutable Field<Type> Ax_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fineBlockAMGLevel(const fineBlockAMGLevel<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const fineBlockAMGLevel<Type>&);


public:

    //- Runtime type information
    TypeName("fineBlockAMGLevel");


    // Constructors

        //- Construct from components
        fineBlockAMGLevel
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const word& coarseningType,
            const label groupSize,
            const label minCoarseEqns
        );


    //- Destructor
    virtual ~fineBlockAMGLevel()
    {}


    // Member Functions

        //- Return reference to dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return reference to x
        virtual Field<Type>& x();

        //- Return reference to b
        virtual Field<Type>& b();

        //- Calculate residual
        virtual void residual
        (
            const Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& res
        ) const;

        //- Restrict residual
        virtual void restrictResidual
        (
            const Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& xBuffer,
            Field<Type>& coarseRes,
            bool preSweepsDone
        ) const;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            Field<Type>& x,
            const Field<Type>& coarseX
        ) const;

        //- Smooth level
        virtual void smooth
        (
            Field<Type>& x,
            const Field<Type>& b,
            const label nSweeps
        ) const;

        //- Solve level
        virtual void solve
        (
            Field<Type>& x,
            const Field<Type>& b,
            const scalar tolerance,
            const scalar relTol
        ) const;

        //- Scale x
        virtual void scaleX
        (
            Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& xBuffer
        ) const;

        //- Create next level from current level
        virtual autoPtr<BlockAMGLevel<Type> > makeNextLevel() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fineBlockAMGLevel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
