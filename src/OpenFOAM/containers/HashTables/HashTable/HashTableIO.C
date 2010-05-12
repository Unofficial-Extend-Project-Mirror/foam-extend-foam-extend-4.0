/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "HashTable.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(Istream& is, const label size)
:
    tableSize_(size),
    table_(new hashedEntry*[tableSize_]),
    nElmts_(0),
    endIter_(*this, NULL, 0),
    endConstIter_(*this, NULL, 0)
{
    for (label i=0; i<tableSize_; i++)
    {
        table_[i] = 0;
    }

    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>(Istream& is, HashTable<T, Key, Hash>& L)
{
    is.fatalCheck("operator>>(Istream&, HashTable<T, Key, Hash>&)");

    // Anull list
    L.clear();

    is.fatalCheck("operator>>(Istream& is, HashTable<T, Key, Hash>& L)");

    token firstToken(is);

    is.fatalCheck
    (
        "operator>>(Istream& is, HashTable<T, Key, Hash>& L) : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char listDelimiter = is.readBeginList("HashTable<T, Key, Hash>");

        if (s)
        {
            if (2*s > L.tableSize_)
            {
                L.resize(2*s);
            }

            if (listDelimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; i++)
                {
                    Key key;
                    is >> key;
                    L.insert(key, pTraits<T>(is));

                    is.fatalCheck
                    (
                        "operator>>(Istream&, HashTable<T, Key, Hash>&) : "
                        "reading entry"
                    );
                }
            }
            else
            {
                FatalIOErrorIn
                (
                    "operator>>(Istream& is, HashTable<T, Key, Hash>& L)",
                    is
                )   << "incorrect first token, '(', found " << firstToken.info()
                    << exit(FatalIOError);
            }
        }

        // Read end of contents
        is.readEndList("HashTable");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "operator>>(Istream& is, HashTable<T, Key, Hash>& L)",
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            Key key;
            is >> key;
            T element;
            is >> element;
            L.insert(key, element);

            is.fatalCheck
            (
                "operator>>(Istream&, HashTable<T, Key, Hash>&) : "
                "reading entry"
            );

            is >> lastToken;
        }
    }
    else
    {
        FatalIOErrorIn
        (
            "operator>>(Istream& is, HashTable<T, Key, Hash>& L)",
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck("operator>>(Istream&, HashTable<T, Key, Hash>&)");

    return is;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::operator<<(Ostream& os, const HashTable<T, Key, Hash>& L)
{
    // Write size of HashTable
    os << nl << L.size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST << nl;

    // Write HashTable contents
    for
    (
        typename HashTable<T, Key, Hash>::const_iterator iter = L.begin();
        iter != L.end();
        ++iter
    )
    {
        os << iter.key() << token::SPACE << iter() << nl;
    }

    // Write end of contents
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const HashTable&)");

    return os;
}


// ************************************************************************* //
