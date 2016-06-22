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

SourceFiles
    tensorBlockILUCpPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockILUCpPrecon_H
#define tensorBlockILUCpPrecon_H

#include "BlockILUCpPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate active type factorization
template<>
template<>
void BlockILUCpPrecon<tensor>::calcActiveTypeFactorization
(
    tensorField& preconD,
    tensorField& extUpper,
    tensorField& extLower
) const;


// Calculate factorization (constructor helper)
template<>
void BlockILUCpPrecon<tensor>::calcFactorization() const;


// Precondition
template<>
void BlockILUCpPrecon<tensor>::precondition
(
    tensorField& x,
    const tensorField& b
) const;


// Precondition transpose
template<>
void BlockILUCpPrecon<tensor>::preconditionT
(
    tensorField& xT,
    const tensorField& bT
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
