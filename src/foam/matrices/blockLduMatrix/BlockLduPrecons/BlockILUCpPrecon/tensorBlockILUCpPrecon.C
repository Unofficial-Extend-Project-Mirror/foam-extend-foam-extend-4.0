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
    BlockILUCpPrecon

Description
    Template specialisation for tensor block ILUCp preconditioning

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockILUCpPrecon_H
#define tensorBlockILUCpPrecon_H

#include "BlockILUCpPrecon.H"
#include "tensorBlockILUCpPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
template<>
void BlockILUCpPrecon<tensor>::calcActiveTypeFactorization
(
    tensorField& preconD,
    tensorField& extUpper,
    tensorField& extLower
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUCpPrecon<tensor>::calcFactorization");
}


template<>
void BlockILUCpPrecon<tensor>::calcFactorization() const
{
    // Get upper and lower matrix factors
    CoeffField<tensor>& Lower = extBlockMatrix_.extendedLower();
    CoeffField<tensor>& Upper = extBlockMatrix_.extendedUpper();

    // Calculate factorization
    // Note: lower, diag and upper must have same type as required by the
    // algorithm. This is handled by lowest possible promotion
    if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asScalar(),
                Upper.asScalar(),
                Lower.asScalar()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asLinear(), // Promotes to linear
                Upper.asLinear(),
                Lower.asLinear()
            );
        }
    }
    else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asLinear(),
                Upper.asLinear(), // Promotes to linear
                Lower.asLinear()  // Promotes to linear
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asLinear(),
                Upper.asLinear(),
                Lower.asLinear()
            );
        }
    }
}


template<>
void BlockILUCpPrecon<tensor>::precondition
(
    tensorField& x,
    const tensorField& b
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUCpPrecon<tensor>::precondition");
}


template<>
void BlockILUCpPrecon<tensor>::preconditionT
(
    tensorField& xT,
    const tensorField& bT
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUCpPrecon<tensor>::preconditionT");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
