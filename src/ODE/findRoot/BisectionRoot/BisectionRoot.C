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
    BisectionRoot

Description
    Bisection root Based on Numerical Recipes in C++, Section 9.1, page 358.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "BisectionRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Func>
const Foam::label Foam::BisectionRoot<Func>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::BisectionRoot<Func>::BisectionRoot(const Func& f, const scalar eps)
:
    f_(f),
    eps_(eps)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Func>
Foam::scalar Foam::BisectionRoot<Func>::root
(
    const scalar x0,
    const scalar x1
) const
{
    scalar f, fMid, dx, rtb, xMid;

    f = f_(x0);
    fMid = f_(x1);

    if (f*fMid >= 0)
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::BisectionRoot<Func>::root\n"
            "(\n"
            "    const scalar x0,\n"
            "    const scalar x1\n"
            ") const"
        )   << "Root is not bracketed.  f(x0) = " << f << " f(x1) = " << fMid
            << abort(FatalError);
    }

    // Orient the search such that f > 0 lies at x + dx
    if (f < 0)
    {
        dx = x1 - x0;
        rtb = x0;
    }
    else
    {
        dx = x0 - x1;
        rtb = x1;
    }

    for (label nIter = 0; nIter < maxIter; nIter++)
    {
        dx *= 0.5;
        xMid = rtb + dx;

        fMid = f_(xMid);

        if (fMid <= 0)
        {
            rtb = xMid;
        }

        if (mag(dx) < eps_ || mag(fMid) < SMALL)
        {
            return rtb;
        }
    }

    FatalErrorIn
    (
        "Foam::scalar Foam::BisectionRoot<Func>::root\n"
        "(\n"
        "    const scalar x0,\n"
        "    const scalar x1\n"
        ") const"
    )   << "Maximum number of iterations exceeded" << abort(FatalError);

    // Dummy return to keep compiler happy
    return x0;
}


// ************************************************************************* //
