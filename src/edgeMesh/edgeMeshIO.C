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

#include "edgeMesh.H"
#include "IFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from file
Foam::edgeMesh::edgeMesh(const fileName& fname)
:
    points_(0),
    edges_(0),
    pointEdgesPtr_(NULL)
{
    IFstream is(fname);

    if (is.good())
    {
        is >> points_ >> edges_;
    }
    else
    {
        FatalErrorIn("edgeMesh::edgeMesh(const fileName&)")
            << "cannot open file " << fname
            << abort(FatalError);
    }
}


// construct from Istream
Foam::edgeMesh::edgeMesh(Istream& is)
:
    points_(is),
    edges_(is),
    pointEdgesPtr_(NULL)
{
    // Check state of Istream
    is.check("edgeMesh::edgeMesh(Istream&)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const edgeMesh& em)
{
    os  << em.points_ << nl << em.edges_ << endl;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const edgeMesh&)");

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, edgeMesh& em)
{
    is >> em.points_ >> em.edges_;

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, edgeMesh&)");

    return is;
}


// ************************************************************************* //
