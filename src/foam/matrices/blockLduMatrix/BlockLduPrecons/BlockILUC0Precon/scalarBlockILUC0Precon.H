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
    BlockILUC0Precon

Description
    Template specialisation for scalar block ILUCp preconditioning

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

SourceFiles
    scalarBlockILUC0Precon.C

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockILUC0Precon_H
#define scalarBlockILUC0Precon_H

#include "BlockILUC0Precon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate active type factorization
template<>
template<>
void BlockILUC0Precon<scalar>::calcActiveTypeFactorization
(
    scalarField& preconD,
    scalarField& preconUpper,
    scalarField& preconLower
) const;


// Calculate factorization (constructor helper)
template<>
void BlockILUC0Precon<scalar>::calcFactorization() const;


// Precondition
template<>
void BlockILUC0Precon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const;


// Precondition transpose
template<>
void BlockILUC0Precon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
