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
    BlockAMGCycle

Description
    Algebraic multigrid fixed cycle class for the BlockLduMatrix

Author
    Klas Jareteg, 2012-12-12

SourceFiles
    BlockAMGCycle.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAMGCycle_H
#define BlockAMGCycle_H

#include "BlockAMGLevel.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class fineBlockAMGLevel;

template<class Type>
class coarseBlockAMGLevel;

template<class Type>
class BlockMatrixCoarsening;

template<class Type>
class BlockAMGLevel;


/*---------------------------------------------------------------------------*\
                      Class blockAMGCycleName Declaration
\*---------------------------------------------------------------------------*/

class blockAMGCycleName
{
public:

    // Public enumerations

        //- Cycle type
        enum cycleType
        {
            V_CYCLE,
            W_CYCLE,
            F_CYCLE
        };

    //- Static data

        //- Cycle type names
        static const NamedEnum<cycleType, 3> cycleNames_;
};


/*---------------------------------------------------------------------------*\
                           Class BlockAMGCycle Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAMGCycle
:
    public blockAMGCycleName
{
    // Private data

        //- Pointer to current AMG level
        autoPtr<BlockAMGLevel<Type> > levelPtr_;

        //- Pointer to coarse AMG cycle
        BlockAMGCycle<Type>* coarseLevelPtr_;

        //- Number of coarse levels
        label nLevels_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockAMGCycle(const BlockAMGCycle<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAMGCycle<Type>&);


public:

    //- Runtime type information
    TypeName("BlockAMGCycle");


    // Constructors

        //- Construct from AMG level
        BlockAMGCycle(autoPtr<BlockAMGLevel<Type> > levelPtr);


    //- Destructor
    virtual ~BlockAMGCycle();


    // Member Functions

        //- Make coarse levels
        void makeCoarseLevels(const label nMaxLevels);

        //- Return number of levels
        label nLevels() const
        {
            return nLevels_;
        }

        //- Calculate residual
        void residual
        (
            const Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& res
        ) const
        {
            levelPtr_->residual(x, b, res);
        }


        //- Fixed cycle
        void fixedCycle
        (
            Field<Type>& x,
            const Field<Type>& b,
            Field<Type>& xBuffer,
            const BlockAMGCycle::cycleType cycle,
            const label nPreSweeps,
            const label nPostSweeps,
            const bool scale
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockAMGCycle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
