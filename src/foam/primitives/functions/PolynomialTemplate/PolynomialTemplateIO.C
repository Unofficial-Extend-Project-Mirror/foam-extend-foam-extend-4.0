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

#include "PolynomialTemplate.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Polynomial<PolySize>& poly
)
{
    os  << static_cast
            <VectorSpace<Polynomial<PolySize>, scalar, PolySize> >(poly);

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Polynomial<PolySize>&)"
    );

    return os;
}


// ************************************************************************* //
