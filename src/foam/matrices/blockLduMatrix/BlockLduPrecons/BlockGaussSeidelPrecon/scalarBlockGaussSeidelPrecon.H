/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    BlockGaussSeidelPrecon

Description
    Template specialisation for scalar block Gauss-Seidel preconditioning

SourceFiles
    scalarBlockGaussSeidelPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockGaussSeidelPrecon_H
#define scalarBlockGaussSeidelPrecon_H

#include "BlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
template<>
void BlockGaussSeidelPrecon<scalar>::BlockSweep
(
    scalarField& x,
    const scalarField& dD,
    const scalarField& upper,
    const scalarField& b
) const;

template<>
template<>
void BlockGaussSeidelPrecon<scalar>::BlockSweep
(
    scalarField& x,
    const scalarField& dD,
    const scalarField& upper,
    const scalarField& lower,
    const scalarField& b
) const;


template<>
void BlockGaussSeidelPrecon<scalar>::precondition
(
    Field<scalar>& x,
    const Field<scalar>& b
) const;


template<>
void BlockGaussSeidelPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
