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
    pamgPolicy

Description
    Pairwise agglomerative AMG policy.  Legacy implementation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    pamgPolicy.C

\*---------------------------------------------------------------------------*/

#ifndef pamgPolicy_H
#define pamgPolicy_H

#include "amgPolicy.H"
#include "lduMatrix.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pamgPolicy Declaration
\*---------------------------------------------------------------------------*/

class pamgPolicy
:
    public amgPolicy
{
    // Private Data

        //- Reference to matrix
        const lduMatrix& matrix_;

        //- Child array: for each fine equation give coarse cluster index
        labelField child_;

        //- Number of solo cells
        label nSolo_;

        //- Number of coarse equations
        label nCoarseEqns_;

        //- Can a coarse level be constructed?
        bool coarsen_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pamgPolicy(const pamgPolicy&);

        //- Disallow default bitwise assignment
        void operator=(const pamgPolicy&);

        //- Calculate child
        void calcChild();


    // Private Static Data

        //- Diagonal scaling factor
        static const debug::tolerancesSwitch diagFactor_;


public:

        //- Runtime type information
        TypeName("PAMG");


    // Constructors

        //- Construct from matrix and group size
        pamgPolicy
        (
            const lduMatrix& matrix,
            const label groupSize,
            const label minCoarseEqns
        );

    // Destructor

        virtual ~pamgPolicy();


    // Member Functions

        //- Can a coarse level be constructed?
        virtual bool coarsen() const
        {
            return coarsen_;
        }

        //- Restrict matrix
        virtual autoPtr<amgMatrix> restrictMatrix
        (
            const FieldField<Field, scalar>& bouCoeffs,
            const FieldField<Field, scalar>& intCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        )  const;

        //- Restrict residual
        virtual void restrictResidual
        (
            const scalarField& res,
            scalarField& coarseRes
        ) const;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            scalarField& x,
            const scalarField& coarseX
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
