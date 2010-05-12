/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-8 H. Jasak All rights reserved
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

        //- Child array
        labelField child_;

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
