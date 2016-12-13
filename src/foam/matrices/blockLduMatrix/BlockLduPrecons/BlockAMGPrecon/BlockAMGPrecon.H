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
    BlockAMGPrecon

Description
    AMG preconditioning

Author
    Klas Jareteg, 2012-12-13

SourceFiles
    BlockAMGPreconI.H
    BlockAMGPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAMGPrecon_H
#define BlockAMGPrecon_H

#include "BlockLduPrecon.H"
#include "Switch.H"
#include "BlockAMGCycle.H"
#include "BlockAMGLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockAMGCycle;

/*---------------------------------------------------------------------------*\
                       Class BlockAMGPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAMGPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Data

        //- Cycle type
        const typename BlockAMGCycle<Type>::cycleType cycle_;

        //- Number of pre-sweeps
        const label nPreSweeps_;

        //- Number of post-sweeps
        const label nPostSweeps_;

        //- Maximum number of levels
        const label nMaxLevels_;

        //- Scaling
        const Switch scale_;

        //- AMG cycle pointer
        autoPtr<BlockAMGCycle<Type> > amgPtr_;

        //- x buffer
        mutable Field<Type> xBuffer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockAMGPrecon(const BlockAMGPrecon<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAMGPrecon<Type>&);


public:

    //- Runtime type information
    TypeName("BlockAMG");


    // Constructors

        //- Construct from matrix and dictionary
        BlockAMGPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BlockAMGPrecon();


    // Member Functions

        //- Return number of AMG levels
        label nLevels() const;

        //- Calculate residual
        const Field<Type>& residual
        (
            const Field<Type>& x,
            const Field<Type>& b
        ) const;

        //- Perform a single cycle
        void cycle
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const;

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockAMGPrecon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
