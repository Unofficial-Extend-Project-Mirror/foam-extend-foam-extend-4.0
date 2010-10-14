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

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockCholeskyPrecon_H
#define tensorBlockCholeskyPrecon_H

#include "BlockCholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockCholeskyPrecon<tensor>::calcPreconDiag()
{
    // Decoupled version
    calcDecoupledPreconDiag();
}


template<>
void Foam::BlockCholeskyPrecon<tensor>::precondition
(
    tensorField& x,
    const tensorField& b
) const
{
    // Decoupled version
    decoupledPrecondition(x, b);
}


template<>
void Foam::BlockCholeskyPrecon<tensor>::preconditionT
(
    tensorField& xT,
    const tensorField& bT
) const
{
    // Decoupled version
    decoupledPreconditionT(xT, bT);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
