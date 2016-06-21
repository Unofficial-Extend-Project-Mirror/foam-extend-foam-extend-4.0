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
    NewtonSecantRoot

Description
    NewtonSecant root finder for better stability.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)

Author
    Aleksandar Jemcov.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "NewtonSecantRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<typename Func, typename Deriv>
const Foam::label Foam::NewtonSecantRoot<Func, Deriv>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Func, typename Deriv>
Foam::NewtonSecantRoot<Func, Deriv>::NewtonSecantRoot
(
    const Func& f,
    const Deriv& d,
    const scalar eps
)
:
    f_(f),
    d_(d),
    eps_(eps)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<typename Func, typename Deriv>
Foam::scalar Foam::NewtonSecantRoot<Func, Deriv>::root
(
    scalar xN
) const
{
    if (0 == d_(xN))
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::NewtonSecantRoot<Func, Deriv>::root\n"
            "(\n"
            "    const scalar xN,\n"
            ") const"
        )   << "Jacobian equal to zero.  f'(xN) = " << d_(xN)
            << abort(FatalError);
    }

    scalar xNp1;

    for (label nIter = 0; nIter < maxIter; ++nIter)
    {
        scalar fN = this->f_(xN);
        scalar dN = this->d_(xN);
        scalar dx = fN/dN;

        scalar xBar = xN - dx;

        xNp1 = xN - fN*fN/(dN*(fN - f_(xBar)));

        if (mag(xN - xNp1) <= this->eps_)
        {
            return xNp1;
        }

        xN = xNp1;
    }

    return xNp1;
}


// ************************************************************************* //
