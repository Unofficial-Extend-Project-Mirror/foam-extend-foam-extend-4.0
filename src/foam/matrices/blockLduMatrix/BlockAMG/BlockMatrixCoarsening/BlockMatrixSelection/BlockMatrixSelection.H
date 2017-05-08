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
    BlockMatrixSelection

Description
    Selective AMG policy for block matrices

Author
    Tessa Uroic, FMENA

SourceFiles
    BlockMatrixSelection.C

\*---------------------------------------------------------------------------*/

#ifndef BlockMatrixSelection_H
#define BlockMatrixSelection_H

#include "BlockMatrixCoarsening.H"
#include "BlockLduMatrix.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeff.H"
#include "crMatrix.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class BlockMatrixSelection Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockMatrixSelection
:
    public BlockMatrixCoarsening<Type>
{
    // Private Data

        //- Reference to matrix
        const BlockLduMatrix<Type>& matrix_;

        //- Reference to a templated norm calculator
        //- INFO: Component norm is recommended, based on the pressure equation
        autoPtr<BlockCoeffNorm<Type> > normPtr_;

        //- Number of coarse equations
        label nCoarseEqns_;

        //- Can a coarse level be constructed?
        bool coarsen_;

        //- Prolongation matrix
        crMatrix* Pptr_;

        //- Restriction matrix
        crMatrix* Rptr_;

        //- Coarsening array - array with labels of coarse/fine equations
        //  (coarse = coarseID / fine = -1)
        labelList rowLabel_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockMatrixSelection(const BlockMatrixSelection<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockMatrixSelection<Type>&);

        //- Calculate restriction and prolongation
        void calcCoarsening();


    // Private enumerations

        // Equation type
        enum equationType
        {
            UNDECIDED = -2,
            FINE = -1,
            COARSE = 0     // Later replaced with coarse matrix equation index
        };


    // Private Static Data

        //- Weighting factor
        static const debug::tolerancesSwitch epsilon_;


public:

        //- Runtime type information
        TypeName("SAMG");


    // Constructors

        //- Construct from matrix and group size
        BlockMatrixSelection
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict,
            const label groupSize,
            const label minCoarseEqns
        );

    //- Destructor
    virtual ~BlockMatrixSelection();


    // Member Functions

        //- Return array with labels of equations
        const labelList& coarseningLabels() const
        {
            return rowLabel_;
        }

        //- Can a coarse level be constructed?
        virtual bool coarsen() const
        {
            return coarsen_;
        }

        //- Restrict matrix
        virtual autoPtr<BlockAMGLevel<Type> > restrictMatrix()  const;

        //- Restrict residual
        virtual void restrictResidual
        (
            const Field<Type>& res,
            Field<Type>& coarseRes
        ) const;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            Field<Type>& x,
            const Field<Type>& coarseX
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockMatrixSelection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
