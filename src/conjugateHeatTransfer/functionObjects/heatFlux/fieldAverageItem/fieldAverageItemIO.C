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

#include "fieldAverageItem.H"
#include "IOstreams.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAverageItem::fieldAverageItem(Istream& is)
:
    fieldName_("unknown"),
    mean_(0),
    prime2Mean_(0)
{
    is.check("Foam::fieldAverageItem::fieldAverageItem(Foam::Istream&)");

    const dictionaryEntry entry(dictionary::null, is);

    fieldName_ = entry.keyword();
    entry.lookup("mean") >> mean_;
    entry.lookup("prime2Mean") >> prime2Mean_;
    base_ = baseTypeNames_[entry.lookup("base")];
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, fieldAverageItem& faItem)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::fieldAverageItem&)"
    );

    const dictionaryEntry entry(dictionary::null, is);

    faItem.fieldName_ = entry.keyword();
    entry.lookup("mean") >> faItem.mean_;
    entry.lookup("prime2Mean") >> faItem.prime2Mean_;
    faItem.base_ = faItem.baseTypeNames_[entry.lookup("base")];

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const fieldAverageItem& faItem)
{
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::fieldAverageItem&)"
    );

    os  << faItem.fieldName_ << nl << token::BEGIN_BLOCK << nl;
    os.writeKeyword("mean") << faItem.mean_ << token::END_STATEMENT << nl;
    os.writeKeyword("prime2Mean") << faItem.mean_
        << token::END_STATEMENT << nl;
    os.writeKeyword("base") << faItem.baseTypeNames_[faItem.base_]
        << token::END_STATEMENT << nl << token::END_BLOCK << nl;

    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::fieldAverageItem&)"
    );

    return os;
}


// ************************************************************************* //
