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
    BlockGaussSeidelSmoother

Description
    Symmetric Gauss-Seidel smoother with prescribed number of smoothing sweeps

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    blockGaussSeidelSmoothers.C

\*---------------------------------------------------------------------------*/

#ifndef BlockGaussSeidelSmoother_H
#define BlockGaussSeidelSmoother_H

#include "BlockLduSmoother.H"
#include "blockGaussSeidelPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class BlockGaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockGaussSeidelSmoother
:
    public BlockLduSmoother<Type>
{
    // Private Data

        //- Gauss-Seidel preconditioner
        BlockGaussSeidelPrecon<Type> gs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockGaussSeidelSmoother(const BlockGaussSeidelSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const BlockGaussSeidelSmoother&);


public:

    //- Runtime type information
    TypeName("GaussSeidel");


    // Constructors

        //- Construct from components
        BlockGaussSeidelSmoother
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        )
        :
            BlockLduSmoother<Type>(matrix),
            gs_(matrix)
        {}


    // Destructor

        virtual ~BlockGaussSeidelSmoother()
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
                gs_.precondition(x, b);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
