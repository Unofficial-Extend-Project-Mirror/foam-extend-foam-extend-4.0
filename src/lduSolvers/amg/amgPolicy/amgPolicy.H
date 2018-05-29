/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

        //- Matrix
        const lduMatrix& matrix_;

        //- Boundary coeffs
        const FieldField<Field, scalar>& bouCoeffs_;

        //- Internal coeffs
        const FieldField<Field, scalar>& intCoeffs_;

        //- Interface fields
        const lduInterfaceFieldPtrsList& interfaceFields_;

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
                const FieldField<Field, scalar>& bouCoeffs,
                const FieldField<Field, scalar>& intCoeffs,
                const lduInterfaceFieldPtrsList& interfaceFields,
                const label groupSize,
                const label minCoarseEqns
            ),
            (
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaceFields,
                groupSize,
                minCoarseEqns
            )
        );


    // Selectors

        //- Select given name, group size and matrix
        static autoPtr<amgPolicy> New
        (
            const word& policyType,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& bouCoeffs,
            const FieldField<Field, scalar>& intCoeffs,
            const lduInterfaceFieldPtrsList& interfaceFields,
            const label groupSize,
            const label minCoarseEqns
        );


    // Constructors

        //- Construct from components
        amgPolicy
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& bouCoeffs,
            const FieldField<Field, scalar>& intCoeffs,
            const lduInterfaceFieldPtrsList& interfaceFields,
            const label groupSize,
            const label minCoarseEqns
        )
        :
            matrix_(matrix),
            bouCoeffs_(bouCoeffs),
            intCoeffs_(intCoeffs),
            interfaceFields_(interfaceFields),
            groupSize_(groupSize),
            minCoarseEqns_(minCoarseEqns)
        {}


    //- Destructor
    virtual ~amgPolicy()
    {}


    // Member Functions

        // Matrix access

            //- Return reference to matrix
            const lduMatrix& matrix() const
            {
                return matrix_;
            }

            //- Return reference to boundary coeffs
            const FieldField<Field, scalar>& bouCoeffs() const
            {
                return bouCoeffs_;
            }

            //- Return reference to internal coeffs
            const FieldField<Field, scalar>& intCoeffs() const
            {
                return intCoeffs_;
            }

            //- Return reference to interface fields
            const lduInterfaceFieldPtrsList& interfaceFields() const
            {
                return interfaceFields_;
            }


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
        virtual autoPtr<amgMatrix> restrictMatrix() const = 0;

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
