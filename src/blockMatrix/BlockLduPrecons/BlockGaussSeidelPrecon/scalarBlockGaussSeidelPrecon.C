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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    BlockGaussSeidelPrecon

Description
    Template specialisation for scalar block Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockGaussSeidelPrecon_H
#define scalarBlockGaussSeidelPrecon_H

#include "BlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockGaussSeidelPrecon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const
{
    if (matrix_.diagonal())
    {
        const scalarField& d = matrix_.diag();

        forAll (x, i)
        {
            x[i] = b[i]/d[i];
        }
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        scalarField dD = 1.0/matrix_.diag();
        const scalarField& LowerCoeff = matrix_.lower();
        const scalarField& UpperCoeff = matrix_.upper();

        BlockSweep
        (
            x,
            dD,
            LowerCoeff,
            UpperCoeff,
            b
        );
    }
}


template<>
void Foam::BlockGaussSeidelPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const
{
    if (matrix_.diagonal())
    {
        const scalarField& d = matrix_.diag();

        forAll (xT, i)
        {
            xT[i] = bT[i]/d[i];
        }
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        scalarField dD = 1.0/matrix_.diag();
        const scalarField& LowerCoeff = matrix_.lower();
        const scalarField& UpperCoeff = matrix_.upper();

        // Swap lower and upper coefficients, transposed matrix
        BlockSweep
        (
            xT,
            dD,
            UpperCoeff,
            LowerCoeff,
            bT
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
