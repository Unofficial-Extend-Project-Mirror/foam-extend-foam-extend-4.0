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

#include "SIBS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::SIBS::polyExtrapolate
(
    const label iest,
    const scalar xest,
    const scalarField& yest,
    scalarField& yz,
    scalarField& dy,
    scalarField& x,
    scalarRectangularMatrix& d
) const
{
    label n = yz.size();

    x[iest] = xest;

    for (register label j=0; j<n; j++)
    {
        dy[j] = yz[j] = yest[j];
    }

    if (iest == 0)
    {
        for (register label j=0; j<n; j++)
        {
            d[j][0] = yest[j];
        }
    }
    else
    {
        scalarField c(yest);

        for (register label k1=0; k1<iest; k1++)
        {
            scalar delta = 1.0/(x[iest - k1 - 1] - xest);
            scalar f1 = xest*delta;
            scalar f2 = x[iest - k1 - 1]*delta;

            for (register label j=0; j<n; j++)
            {
                scalar q = d[j][k1];
                d[j][k1] = dy[j];
                delta = c[j] - q;
                dy[j] = f1*delta;
                c[j] = f2*delta;
                yz[j] += dy[j];
            }
        }

        for (register label j=0; j<n; j++)
        {
            d[j][iest] = dy[j];
        }
    }
}


// ************************************************************************* //
