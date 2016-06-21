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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::STLtriangle::STLtriangle()
{}


inline Foam::STLtriangle::STLtriangle
(
    const STLpoint& normal,
    const STLpoint& a,
    const STLpoint& b,
    const STLpoint& c,
    unsigned short attrib
)
:
    normal_(normal),
    a_(a),
    b_(b),
    c_(c),
    attrib_(attrib)
{}


inline Foam::STLtriangle::STLtriangle(istream& is)
{
    read(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::STLpoint& Foam::STLtriangle::normal() const
{
    return normal_;
}


inline const Foam::STLpoint& Foam::STLtriangle::a() const
{
    return a_;
}


inline const Foam::STLpoint& Foam::STLtriangle::b() const
{
    return b_;
}


inline const Foam::STLpoint& Foam::STLtriangle::c() const
{
    return c_;
}


inline unsigned short Foam::STLtriangle::attrib() const
{
    return attrib_;
}


inline void Foam::STLtriangle::read(istream& is)
{
    is.read(reinterpret_cast<char*>(this), 4*sizeof(STLpoint));
    is.read(reinterpret_cast<char*>(&attrib_), sizeof(STLattrib));
}


inline void Foam::STLtriangle::write(ostream& os)
{
    os.write(reinterpret_cast<char*>(this), 4*sizeof(STLpoint));
    os.write(reinterpret_cast<char*>(&attrib_), sizeof(STLattrib));
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const STLtriangle& t)
{
    os  << t.normal_ << token::SPACE
        << t.a_ << token::SPACE
        << t.b_ << token::SPACE
        << t.c_ << token::SPACE
        << t.attrib_;

    return os;
}


// ************************************************************************* //
