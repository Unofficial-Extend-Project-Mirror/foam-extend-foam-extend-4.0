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

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"
#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance
)
{
    if (blockLduMatrix::debug >= 2)
    {
        Info<< solverName_
            << ":  Iteration " << nIterations_
            << " residual = " << finalResidual_
            << " mag = " << mag(finalResidual_)
            << " tol = "
            << Foam::max(Tolerance, RelTolerance*mag(initialResidual_))
            << endl;
    }

    // Reconsider evaluation of the final residual residualVector
    // - mag(residualVector) = sqrt(sum(sqr(cmpt))).  Currently used - strict
    // - cmptSum(residualVector) = consistent with 1-norm
    // - cmptMax(residualVector) = consistent with inftyNorm
    // HJ, 29/May/2014
    if
    (
        mag(finalResidual_) < Tolerance
     || (
            RelTolerance > SMALL
         && mag(finalResidual_) <= RelTolerance*mag(initialResidual_)
        )
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkSingularity
(
    const scalar& residual
)
{
    if (mag(residual) > VSMALL)
    {
        singular_ = false;
    }
    else
    {
        singular_ = true;
    }

    return singular_;
}


template<class Type>
void Foam::BlockSolverPerformance<Type>::print() const
{
    if (blockLduMatrix::debug)
    {
        Info<< solverName_ << ":  Solving for " << fieldName_;

        if (singular())
        {
            Info<< ":  solution singularity" << endl;
        }
        else
        {
            Info<< ", Initial residual = " << initialResidual_
                << ", Final residual = " << finalResidual_
                << ", No Iterations " << nIterations_
                << endl;
        }
    }
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::operator!=
(
    const BlockSolverPerformance<Type>& bsp
) const
{
    return
    (
        solverName()      != bsp.solverName()
     || fieldName()       != bsp.fieldName()
     || initialResidual() != bsp.initialResidual()
     || finalResidual()   != bsp.finalResidual()
     || nIterations()     != bsp.nIterations()
     || converged()       != bsp.converged()
     || singular()        != bsp.singular()
    );
}


template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    typename Foam::BlockSolverPerformance<Type>& bsp
)
{
    is.readBeginList("BlockSolverPerformance<Type>");
    is  >> bsp.solverName_
        >> bsp.fieldName_
        >> bsp.initialResidual_
        >> bsp.finalResidual_
        >> bsp.nIterations_
        >> bsp.converged_
        >> bsp.singular_;
    is.readEndList("BlockSolverPerformance<Type>");

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const typename Foam::BlockSolverPerformance<Type>& bsp
)
{
    os  << token::BEGIN_LIST
        << bsp.solverName_ << token::SPACE
        << bsp.fieldName_ << token::SPACE
        << bsp.initialResidual_ << token::SPACE
        << bsp.finalResidual_ << token::SPACE
        << bsp.nIterations_ << token::SPACE
        << bsp.converged_ << token::SPACE
        << bsp.singular_ << token::SPACE
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
