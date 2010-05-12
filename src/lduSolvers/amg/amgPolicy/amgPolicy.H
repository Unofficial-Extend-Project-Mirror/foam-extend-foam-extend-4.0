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
    amgPolicy

Description
    Virtual base class for AMG coarsening policy

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    amgPolicy.C

\*---------------------------------------------------------------------------*/

#ifndef amgPolicy_H
#define amgPolicy_H

#include "runTimeSelectionTables.H"
#include "primitiveFieldsFwd.H"
#include "FieldFields.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMatrix;
class amgMatrix;

/*---------------------------------------------------------------------------*\
                           Class amgPolicy Declaration
\*---------------------------------------------------------------------------*/

class amgPolicy
{
    // Private Data

        //- Group size
        label groupSize_;

        //- Minimum number of coarse level equations
        label minCoarseEqns_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        amgPolicy(const amgPolicy&);

        //- Disallow default bitwise assignment
        void operator=(const amgPolicy&);


public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                amgPolicy,
                matrix,
                (
                    const lduMatrix& matrix,
                    const label groupSize,
                    const label minCoarseEqns
                ),
                (matrix, groupSize, minCoarseEqns)
            );


    // Selectors

        //- Select given name, group size and matrix
        static autoPtr<amgPolicy> New
        (
            const word& policyType,
            const lduMatrix& matrix,
            const label groupSize,
            const label minCoarseEqns
        );


    // Constructors

        //- Construct from components
        amgPolicy
        (
            const label groupSize,
            const label minCoarseEqns
        )
        :
            groupSize_(groupSize),
            minCoarseEqns_(minCoarseEqns)
        {}


    // Destructor

        virtual ~amgPolicy()
        {}


    // Member Functions

        //- Return group size
        label groupSize() const
        {
            return groupSize_;
        }

        //- Return minimum number of coarse level equations
        label minCoarseEqns() const
        {
            return minCoarseEqns_;
        }

        //- Can a coarse level be constructed?
        virtual bool coarsen() const = 0;

        //- Restrict matrix
        virtual autoPtr<amgMatrix> restrictMatrix
        (
            const FieldField<Field, scalar>& bouCoeffs,
            const FieldField<Field, scalar>& intCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        ) const = 0;

        //- Restrict residual
        virtual void restrictResidual
        (
            const scalarField& res,
            scalarField& coarseRes
        ) const = 0;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            scalarField& x,
            const scalarField& coarseX
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
