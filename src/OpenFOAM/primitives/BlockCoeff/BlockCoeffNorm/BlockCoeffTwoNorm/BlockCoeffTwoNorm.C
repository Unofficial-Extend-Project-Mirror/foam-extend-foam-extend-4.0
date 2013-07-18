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
    BlockCoeffTwoNorm

\*---------------------------------------------------------------------------*/

#include "BlockCoeffTwoNorm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeffTwoNorm<Type>::BlockCoeffTwoNorm
(
    const dictionary& dict
)
:
    BlockCoeffNorm<Type>(dict),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::BlockCoeffTwoNorm<Type>::normalize
(
    const Foam::BlockCoeff<Type>& a
)
{
    if (a.activeType() == Foam::BlockCoeff<Type>::SCALAR)
    {
        return mag(a.asScalar());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::LINEAR)
    {
        return mag(a.asLinear());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::SQUARE)
    {
        return mag(a.asSquare());
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffTwoNorm<Type>(const BlockCoeff<Type>& a)"
        )   << "Unknown type" << abort(FatalError);

        return 0;
    }

    // Dummy return
    return scalar(0);
}


template<class Type>
void Foam::BlockCoeffTwoNorm<Type>::coeffMag
(
    const Foam::CoeffField<Type>& a,
    Foam::Field<scalar>& b
)
{
    if (a.activeType() == Foam::BlockCoeff<Type>::SCALAR)
    {
        b = mag(a.asScalar());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::LINEAR)
    {
        b = mag(a.asLinear());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::SQUARE)
    {
        b = mag(a.asSquare());
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffTwoNorm<Type>(const BlockCoeff<Type>& b)"
        )   << "Unknown type" << abort(FatalError);
    }
}


// ************************************************************************* //
