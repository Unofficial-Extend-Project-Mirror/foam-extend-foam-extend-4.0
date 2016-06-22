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

#include "polynomial.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polynomial, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polynomial::polynomial(const word& entryName, Istream& is)
:
    DataEntry<scalar>(entryName),
    coeffs_(is)
{
    if (!coeffs_.size())
    {
        FatalErrorIn("Foam::polynomial::polynomial(const word&, Istream&)")
            << "polynomial coefficients for entry " << this->name_
            << " is invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::polynomial::polynomial(const polynomial& poly)
:
    DataEntry<scalar>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polynomial::~polynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::polynomial::value(const scalar x) const
{
    scalar y = 0.0;
    forAll(coeffs_, i)
    {
        y += coeffs_[i].first()*pow(x, coeffs_[i].second());
    }

    return y;
}


Foam::scalar Foam::polynomial::integrate(const scalar x1, const scalar x2) const
{
    scalar intx = 0.0;

    forAll(coeffs_, i)
    {
        intx +=
            coeffs_[i].first()/(coeffs_[i].second() + 1)
           *(
                pow(x2, coeffs_[i].second() + 1)
              - pow(x1, coeffs_[i].second() + 1)
            );
    }

    return intx;
}


// ************************************************************************* //
