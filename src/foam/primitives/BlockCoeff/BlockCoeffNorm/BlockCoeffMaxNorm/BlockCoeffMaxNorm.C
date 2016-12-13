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
    BlockCoeffMaxNorm

\*---------------------------------------------------------------------------*/

#include "BlockCoeffMaxNorm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeffMaxNorm<Type>::BlockCoeffMaxNorm
(
    const dictionary& dict
)
:
    BlockCoeffNorm<Type>(dict),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::BlockCoeffMaxNorm<Type>::normalize
(
    const Foam::BlockCoeff<Type>& a
)
{
    if (a.activeType() == BlockCoeff<Type>::SCALAR)
    {
        return mag(a.asScalar());
    }
    else if (a.activeType() == BlockCoeff<Type>::LINEAR)
    {
        return cmptMax(cmptMag(a.asLinear()));
    }
    else if (a.activeType() == BlockCoeff<Type>::SQUARE)
    {
        return cmptMax(cmptMag(a.asSquare()));
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffMaxNorm<Type>(const BlockCoeff<Type>& b)"
        )   << "Unknown type" << abort(FatalError);

        return 0;
    }
}


template<class Type>
void Foam::BlockCoeffMaxNorm<Type>::coeffMag
(
    const Foam::CoeffField<Type>& a,
    Foam::Field<scalar>& b
)
{
    if (a.activeType() == BlockCoeff<Type>::SCALAR)
    {
        b = mag(a.asScalar());
    }
    else if (a.activeType() == BlockCoeff<Type>::LINEAR)
    {
        b = cmptMax(cmptMag(a.asLinear()));
    }
    else if (a.activeType() == BlockCoeff<Type>::SQUARE)
    {
        b = cmptMax(cmptMag(a.asSquare()));
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffMaxNorm<Type>(const BlockCoeff<Type>& b)"
        )   << "Unknown type" << abort(FatalError);
    }
}


// ************************************************************************* //
