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

#include "geomDecomp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geomDecomp::geomDecomp
(
    const dictionary& decompositionDict,
    const word& derivedType
)
:
    decompositionMethod(decompositionDict),
    geomDecomDict_(decompositionDict.subDict(derivedType + "Coeffs")),
    n_(geomDecomDict_.lookup("n")),
    delta_(readScalar(geomDecomDict_.lookup("delta"))),
    rotDelta_(I)
{
    // Check that the case makes sense
    if (nProcessors_ != n_.x()*n_.y()*n_.z())
    {
        FatalErrorIn
        (
            "geomDecomp::geomDecomp"
            "(const dictionary& decompositionDict)"
        )   << "Wrong number of processor divisions in geomDecomp:" << nl
            << "Number of domains    : " << nProcessors_ << nl
            << "Wanted decomposition : " << n_
            << exit(FatalError);
    }

    scalar d = 1 - 0.5*delta_*delta_;
    scalar d2 = sqr(d);

    scalar a = delta_;
    scalar a2 = sqr(a);

    rotDelta_ = tensor
    (
        d2,         -a*d,         a,
        a*d - a2*d,  a*a2 + d2,  -2*a*d,
        a*d2 + a2,   a*d - a2*d,  d2 - a2
    );
}


// ************************************************************************* //
