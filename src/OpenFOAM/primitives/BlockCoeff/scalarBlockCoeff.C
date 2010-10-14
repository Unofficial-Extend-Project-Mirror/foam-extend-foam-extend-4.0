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
