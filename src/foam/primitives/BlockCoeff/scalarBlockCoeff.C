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

\*---------------------------------------------------------------------------*/

#include "scalarBlockCoeff.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BlockCoeff<Foam::scalar>::BlockCoeff()
:
    scalarCoeff_(pTraits<scalar>::zero)
{}


Foam::BlockCoeff<Foam::scalar>::BlockCoeff(const BlockCoeff<scalar>& f)
:
    scalarCoeff_(f.scalarCoeff_)
{}


Foam::BlockCoeff<Foam::scalar>::BlockCoeff(Istream& is)
:
    scalarCoeff_(readScalar(is))
{}


Foam::BlockCoeff<Foam::scalar> Foam::BlockCoeff<Foam::scalar>::clone() const
{
    return BlockCoeff<scalar>(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BlockCoeff<Foam::scalar>::~BlockCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockCoeffBase::activeLevel
Foam::BlockCoeff<Foam::scalar>::activeType() const
{
    return blockCoeffBase::SCALAR;
}


Foam::scalar Foam::BlockCoeff<Foam::scalar>::component(const direction) const
{
    return scalarCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::BlockCoeff<Foam::scalar>::operator=(const BlockCoeff<scalar>& f)
{
    if (this == &f)
    {
        FatalErrorIn
        (
            "BlockCoeff<scalar>::operator=(const BlockCoeff<scalar>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    scalarCoeff_ = f.scalarCoeff_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const BlockCoeff<scalar>& f)
{
    os << f.scalarCoeff_;

    return os;
}


// ************************************************************************* //
