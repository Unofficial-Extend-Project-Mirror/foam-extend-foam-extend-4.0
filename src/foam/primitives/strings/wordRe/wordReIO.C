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

#include "wordRe.H"
#include "IOstreams.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::wordRe::wordRe(Istream& is)
:
    word(),
    re_(NULL)
{
    is >> *this;
}


Foam::Istream& Foam::operator>>(Istream& is, wordRe& w)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isWord())
    {
        w = t.wordToken();
    }
    else if (t.isString())
    {
        // Auto-tests for regular expression
        w = t.stringToken();
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, wordRe&)", is)
            << "wrong token type - expected word or string found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, wordRe&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const wordRe& w)
{
    os.writeQuoted(w, w.isPattern());
    os.check("Ostream& operator<<(Ostream&, const wordRe&)");
    return os;
}


Foam::Ostream& Foam::wordRe::info(Ostream& os) const
{
    if (isPattern())
    {
        os  << "wordRe(regex) " << *this;
    }
    else
    {
        os  << "wordRe(plain) \"" << *this << '"';
    }
    os.flush();

    return os;
}


// ************************************************************************* //
