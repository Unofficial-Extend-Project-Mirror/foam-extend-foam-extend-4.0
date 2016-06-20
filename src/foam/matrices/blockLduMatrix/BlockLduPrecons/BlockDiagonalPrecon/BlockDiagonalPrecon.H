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
    BlockDiagonalPrecon

Description
    Diagonal preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    blockDiagonalPrecons.C
    scalarBlockDiagonalPrecon.C
    tensorBlockDiagonalPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockDiagonalPrecon_H
#define BlockDiagonalPrecon_H

#include "BlockLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockDiagonalPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockDiagonalPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockDiagonalPrecon(const BlockDiagonalPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockDiagonalPrecon&);


public:

    //- Runtime type information
    TypeName("diagonal");


    // Constructors

        //- Construct from components
        BlockDiagonalPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        )
        :
            BlockLduPrecon<Type>(matrix)
        {}


    // Destructor

        virtual ~BlockDiagonalPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const;

        //- Execute preconditioning on a transposed matrix
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const
        {
            return precondition(xT, bT);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockDiagonalPrecon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
