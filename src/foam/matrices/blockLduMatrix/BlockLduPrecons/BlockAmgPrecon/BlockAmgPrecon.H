/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    BlockAmgPrecon

Description
    AMG preconditioning

Author
    Klas Jareteg, 2012-12-13

SourceFiles
    BlockAmgPreconI.H
    BlockAmgPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAmgPrecon_H
#define BlockAmgPrecon_H

#include "BlockLduPrecon.H"
#include "Switch.H"
#include "BlockAmgCycle.H"
#include "BlockAmgLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockAmgCycle;


/*---------------------------------------------------------------------------*\
                           Class BlockAmgPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAmgPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Data

        //- Cycle type
        const typename BlockAmgCycle<Type>::cycleType cycle_;

        //- Number of pre-sweeps
        const label nPreSweeps_;

        //- Number of post-sweeps
        const label nPostSweeps_;

        //- Maximum number of levels
        const label nMaxLevels_;

        //- Scaling
        const Switch scale_;

        //- AMG cycle pointer
        autoPtr<BlockAmgCycle<Type> > amgPtr_;

        //- x buffer
        mutable Field<Type> xBuffer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockAmgPrecon(const BlockAmgPrecon<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAmgPrecon<Type>&);


public:

    //- Runtime type information
    TypeName("BlockAMG");


    // Constructors

        //- Construct from matrix and dictionary
        BlockAmgPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockAmgPrecon();


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
#   include "BlockAmgPrecon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
