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

Description
    Update of block interfaces

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockLduMatrix<Type>::initInterfaces
(
    TypeField& Ax,
    const TypeField& x
) const
{
    forAll (interfaces_, interfaceI)
    {
        interfaces_[interfaceI].initMatrixUpdate(*this, Ax, x);
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::updateInterfaces
(
    const FieldField<CoeffField, Type>& coeffs,
    TypeField& Ax,
    const TypeField& x
) const
{
    forAll (interfaces_, interfaceI)
    {
        interfaces_[interfaceI].updateMatrix(*this, coeffs, Ax, x);
    }

}


// ************************************************************************* //
