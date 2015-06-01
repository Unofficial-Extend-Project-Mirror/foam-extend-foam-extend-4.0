/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline Foam::SRIFallOffFunction::SRIFallOffFunction
(
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d,
    const scalar e
)
:
    a_(a),
    b_(b),
    c_(c),
    d_(d),
    e_(e)
{}


// Construct from Istream
inline Foam::SRIFallOffFunction::SRIFallOffFunction(Istream& is)
:
    a_(readScalar(is.readBegin("SRIFallOffFunction(Istream&)"))),
    b_(readScalar(is)),
    c_(readScalar(is)),
    d_(readScalar(is)),
    e_(readScalar(is))
{
    is.readEnd("SRIFallOffFunction(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::SRIFallOffFunction::operator()
(
    const scalar T,
    const scalar Pr
) const
{
    scalar X = 1.0/(1.0 + sqr(log10(max(Pr, SMALL))));
    return d_*pow(a_*exp(-b_/T) + exp(-T/c_), X)*pow(T, e_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::SRIFallOffFunction& srifof
)
{
    os  << token::BEGIN_LIST
        << srifof.a_
        << token::SPACE << srifof.b_
        << token::SPACE << srifof.c_
        << token::SPACE << srifof.d_
        << token::SPACE << srifof.e_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
