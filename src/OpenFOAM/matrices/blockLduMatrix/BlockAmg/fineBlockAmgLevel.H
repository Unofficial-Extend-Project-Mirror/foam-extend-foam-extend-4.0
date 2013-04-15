/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    fineBlockAmgLevel

Description
    Fine AMG level refers to matrix, x and b data from the solver or
    preconditioner, implemented for the BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

SourceFiles
    fineBlockAmgLevel.C

\*---------------------------------------------------------------------------*/

#ifndef fineBlockAmgLevel_H
#define fineBlockAmgLevel_H

#include "BlockAmgLevel.H"
#include "BlockLduSmoother.H"
#include "BlockLduMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class coarseBlockAmgLevel;

template<class Type>
class BlockAmgPolicy;

/*---------------------------------------------------------------------------*\
                        Class fineBlockAmgLevel Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fineBlockAmgLevel
:
    public BlockAmgLevel<Type>
{
    // Private data

         //- Reference to matrix
        const BlockLduMatrix<Type>& matrix_;

        //- Dictionary
        const dictionary& dict_;

        //- AMG policy
        autoPtr<BlockAmgPolicy<Type> > policyPtr_;

        //- Smoother
        autoPtr<BlockLduSmoother<Type> > smootherPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fineBlockAmgLevel(const fineBlockAmgLevel<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const fineBlockAmgLevel<Type>&);


public:
    
    //- Runtime type information
        TypeName("fineBlockAmgLevel");


        //- Construct from components
        fineBlockAmgLevel
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const word& policyType,
            const label groupSize,
            const label minCoarseEqns,
            const word& smootherType
        );


    // Destructor

        virtual ~fineBlockAmgLevel()
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
        virtual autoPtr<BlockAmgLevel<Type> > makeNextLevel() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fineBlockAmgLevel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
