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
    BlockCoeffComponentNorm

\*---------------------------------------------------------------------------*/

#include "BlockCoeffComponentNorm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeffComponentNorm<Type>::BlockCoeffComponentNorm
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
Foam::scalar Foam::BlockCoeffComponentNorm<Type>::normalize
(
    const Foam::BlockCoeff<Type>& a
)
{
    return mag(a.component(cmpt_));
}


template<class Type>
void Foam::BlockCoeffComponentNorm<Type>::coeffMag
(
    const Foam::CoeffField<Type>& a,
    Foam::Field<scalar>& b
)
{
    b = mag(a.component(cmpt_));
}


// ************************************************************************* //
