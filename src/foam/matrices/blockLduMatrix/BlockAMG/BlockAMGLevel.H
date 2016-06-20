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
    BlockAMGLevel

Description
    Algebraic multigrid level virtual base class for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

SourceFiles
    BlockAMGLevel.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAMGLevel_H
#define BlockAMGLevel_H

#include "BlockLduMatrix.H"
#include "BlockMatrixCoarsening.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class fineBlockAMGLevel;

template<class Type>
class coarseBlockAMGLevel;

template<class Type>
class BlockAMGCycle;

template<class Type>
class BlockMatrixCoarsening;

template<class Type>
class BlockAMGSolver;

template<class Type>
class BlockAMGPrecon;

/*---------------------------------------------------------------------------*\
                           Class BlockAMGLevel Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAMGLevel
{
public:

    // Runtime type information
    TypeName("BlockAMGLevel");


    //- Destructor
    virtual ~BlockAMGLevel()
    {}


    // Member Functions

        //- Return reference to x
        virtual Field<Type>& x() = 0;

        //- Return reference to b
        virtual Field<Type>& b() = 0;

        //- Calculate residual
        virtual void residual
        (
            const Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& res
        ) const = 0;

        //- Restrict residual
        virtual void restrictResidual
        (
            const Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& xBuffer,
            Field<Type>& coarseRes,
            bool preSweepsDone
        ) const = 0;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            Field<Type>& x,
            const Field<Type>& coarseX
        ) const = 0;

        //- Smooth level
        virtual void smooth
        (
            Field<Type>& x,
            const Field<Type>& b,
            const label nSweeps
        ) const = 0;

        //- Solve level
        virtual void solve
        (
            Field<Type>& x,
            const Field<Type>& b,
            const scalar tolerance,
            const scalar relTol
        ) const = 0;

        //- Scale x
        virtual void scaleX
        (
            Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& xBuffer
        ) const = 0;

        //- Create next level from current level
        virtual autoPtr<BlockAMGLevel> makeNextLevel() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
