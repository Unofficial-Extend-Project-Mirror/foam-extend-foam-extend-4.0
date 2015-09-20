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

#include "equation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Friend IOstream Operators * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, equation& I)
{
    // Acceptable istream formats:
    //      scalar;
    //      string;
    //      [dimensionSet] scalar;
    //      [dimensionSet] string;
    //      name [dimensionSet] scalar;
    //      name [dimensionSet] string;

    if (I.equationName_ == word::null)
    {
        I.equationName_ = "fromIstream";
    }
    token t(is);
    if (t.isString())
    {
        I.rawText_ = t.stringToken();
    }
    else if (t.isNumber())
    {
        I.rawText_ = string(name(t.number()));
    }
    else if (t.isPunctuation())
    {
        is.putBack(t);
        I.changeDimensions_ = true;
        I.overrideDimensions_.reset(dimensionSet(is));
        token t2(is);
        if (t2.isString())
        {
            is.putBack(t2);
            I.rawText_ = string(is);
        }
        else // number
        {
            I.rawText_ = string(name(t.number()));
        }
    }
    else if (t.isWord())
    {
        word garbage(t.wordToken());
        I.changeDimensions_ = true;
        I.overrideDimensions_.reset(dimensionSet(is));
        token t2(is);
        if (t2.isString())
        {
            is.putBack(t2);
            I.rawText_ = string(is);
        }
        else // number
        {
            I.rawText_ = string(name(t.number()));
        }
    }
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const equation& I)
{
    if (I.changeDimensions_)
    {
        os  << I.overrideDimensions_ << token::TAB
            << I.rawText_;
    }
    else
    {
        os << I.rawText_;
    }
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
