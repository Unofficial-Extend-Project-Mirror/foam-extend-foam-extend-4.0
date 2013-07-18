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
    BlockAmgCycle

Description
    Algebraic multigrid fixed cycle class for the BlockLduMatrix

Author
    Klas Jareteg, 2012-12-12

SourceFiles
    BlockAmgCycle.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAmgCycle_H
#define BlockAmgCycle_H

#include "BlockAmgLevel.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class fineBlockAmgLevel;

template<class Type>
class coarseBlockAmgLevel;

template<class Type>
class BlockAmgPolicy;

template<class Type>
class BlockAmgLevel;


/*---------------------------------------------------------------------------*\
                   Class blockAmgCycleName Declaration
\*---------------------------------------------------------------------------*/

class blockAmgCycleName
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
                           Class BlockAmgCycle Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAmgCycle
:
    public blockAmgCycleName
{
    // Private data

        //- Pointer to current AMG level
        autoPtr<BlockAmgLevel<Type> > levelPtr_;

        //- Pointer to coarse AMG cycle
        BlockAmgCycle<Type>* coarseLevelPtr_;

        //- Number of coarse levels
        label nLevels_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockAmgCycle(const BlockAmgCycle<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAmgCycle<Type>&);


public:

    //- Runtime type information
    TypeName("BlockAmgCycle");


    // Constructors

        //- Construct from AMG level
        BlockAmgCycle(autoPtr<BlockAmgLevel<Type> > levelPtr);


    // Destructor

        virtual ~BlockAmgCycle();


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
            const BlockAmgCycle::cycleType cycle,
            const label nPreSweeps,
            const label nPostSweeps,
            const bool scale
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockAmgCycle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
