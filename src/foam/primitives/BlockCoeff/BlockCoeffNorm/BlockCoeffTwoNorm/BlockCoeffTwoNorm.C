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
    dict_(dict),
    cmpt_(readInt(this->dict().lookup("normComponent")))
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
        // Note: for two-norm, use mag
        return sign(a.component(cmpt_))*mag(a.asScalar());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::LINEAR)
    {
        // Note: for two-norm, use mag
        return sign(a.component(cmpt_))*mag(a.asLinear());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::SQUARE)
    {
        // Note: for two-norm, use mag of diag
        typedef typename BlockCoeff<Type>::linearType linearType;

        linearType diag;
        contractLinear(diag, a.asSquare());

        return sign(diag.component(cmpt_))*mag(diag);
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffTwoNorm<Type>::normalize"
            "(const BlockCoeff<Type>& a)"
        )   << "Unknown type" << abort(FatalError);

        return 0;
    }

    // Dummy return
    return 0;
}


template<class Type>
void Foam::BlockCoeffTwoNorm<Type>::normalize
(
    Field<scalar>& b,
    const CoeffField<Type>& a
)
{
    if (a.activeType() == Foam::BlockCoeff<Type>::SCALAR)
    {
        // Note: for two-norm, use mag
        // However, the norm is signed
        // b = mag(a.asScalar());
        b = a.asScalar();
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::LINEAR)
    {
        // Note: for two-norm, use mag = sqrt(sum magSqr(cmpt))
        // However, the norm is signed
        b = sign(a.component(cmpt_))*mag(a.asLinear());
    }
    else if (a.activeType() == Foam::BlockCoeff<Type>::SQUARE)
    {
        // Note: for two-norm, use mag of diag
        typedef typename BlockCoeff<Type>::linearTypeField linearTypeField;

        linearTypeField diag(a.size());
        contractLinear(diag, a.asSquare());

        b = sign(diag.component(cmpt_))*mag(diag);
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
