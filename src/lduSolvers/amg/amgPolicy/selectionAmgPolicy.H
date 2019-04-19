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
    selectionAmgPolicy

Description
    Classical AMG coarsening algorithm. The choice of coarse
    equations is based on the number of strong connections to other
    equations. For optimal performance, strong connections should have sign
    opposite to diagonal coefficient. Coarse level matrix is obtained by
    multiplying the restriction matrix, fine level matrix and prolongation
    matrix. The algorithm closely follows theoretical background from work of
    Klaus Stueben.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

References
    [1] U. Trottenberg, C. Oosterlee, A. Schueller:
            "Multigrid", Elsevier, Academic Press, 2001

SourceFiles
    selectionAmgPolicy.C

\*---------------------------------------------------------------------------*/

#ifndef selectionAmgPolicy_H
#define selectionAmgPolicy_H

#include "amgPolicy.H"
#include "lduMatrix.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class crMatrix;

/*---------------------------------------------------------------------------*\
                      Class selectionAmgPolicy Declaration
\*---------------------------------------------------------------------------*/

class selectionAmgPolicy
:
    public amgPolicy
{
    // Private Data

        //- Number of coarse equations
        label nCoarseEqns_;

        //- Can a coarse level be constructed?
        bool coarsen_;

        //- Prolongation matrix
        crMatrix* Pptr_;

        //- Restriction matrix
        crMatrix* Rptr_;

        //- Coarsening array - array with labels of coarse/fine equations
        //  (coarse = coarseID / fine = FINE (-1)).  See enumeration below
        labelList rowLabel_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        selectionAmgPolicy(const selectionAmgPolicy&);

        //- Disallow default bitwise assignment
        void operator=(const selectionAmgPolicy&);

        //- Filter the prolongation matrix and return only the rows
        //  corresponding to interface
        //- Take the big prolongation matrix and faceCells as arguments
        autoPtr<crMatrix> filterProlongation
        (
            const crMatrix& prolongationMatrix,
            const labelList& fineFaceCells
        ) const;

        //- Calculate restriction and prolongation
        void calcCoarsening();


    // Private enumerations

        //- Equation selection enumeration
        enum
        {
            UNDECIDED = -2,
            FINE = -1,
            COARSE = 0     // Later replaced with coarse matrix equation index
        };


    // Private Static Data

        //- Weighting factor
        static const debug::tolerancesSwitch epsilon_;

        //- Diagonal scaling factor
        static const debug::tolerancesSwitch diagFactor_;


public:

        //- Runtime type information
        TypeName("selectionAMG");


    // Constructors

        //- Construct from matrix and group size
        selectionAmgPolicy
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& bouCoeffs,
            const FieldField<Field, scalar>& intCoeffs,
            const lduInterfaceFieldPtrsList& interfaceFields,
            const label groupSize,
            const label minCoarseEqns
        );


    //- Destructor
    virtual ~selectionAmgPolicy();


    // Member Functions

        //- Can a coarse level be constructed?
        virtual bool coarsen() const
        {
            return coarsen_;
        }

        //- Restrict matrix
        virtual autoPtr<amgMatrix> restrictMatrix()  const;

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

        //- Return array with labels of FINE/COARSE cells
        // TU, 17/May/2017, debugging
        virtual labelList coarseningLabels() const
        {
            return rowLabel_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
