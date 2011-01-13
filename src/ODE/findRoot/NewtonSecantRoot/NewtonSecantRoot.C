/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
const Foam::label Foam::NewtonSecantRoot<Func,Deriv>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Func, typename Deriv>
Foam::NewtonSecantRoot<Func,Deriv>::NewtonSecantRoot
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
Foam::scalar Foam::NewtonSecantRoot<Func,Deriv>::root
(
    scalar xN
) const
{
    if (0 == d_(xN))
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::NewtonSecantRoot<Func,Deriv>::root\n"
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
