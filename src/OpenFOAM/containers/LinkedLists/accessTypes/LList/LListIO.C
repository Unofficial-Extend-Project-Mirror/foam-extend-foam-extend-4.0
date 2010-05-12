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

Description

\*---------------------------------------------------------------------------*/

#include "LList.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Istream
template<class LListBase, class T>
LList<LListBase, T>::LList(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Istream& operator>>(Istream& is, LList<LListBase, T>& L)
{
    // Anull list
    L.clear();

    is.fatalCheck(" operator>>(Istream& is, LList<LListBase, T>& L)");

    token firstToken(is);

    is.fatalCheck
    (
        " operator>>(Istream& is, LList<LListBase, T>& L) : reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char listDelimiter = is.readBeginList("LList<LListBase, T>");

        if (s)
        {
            if (listDelimiter == token::BEGIN_LIST)
            {
                for (register label i=0; i<s; i++)
                {
                    T element;
                    is >> element;
                    L.append(element);
                }
            }
            else
            {
                T element;
                is >> element;

                for (register label i=0; i<s; i++)
                {
                    L.append(element);
                }
            }
        }

        // Read end of contents
        is.readEndList("LList");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                " operator>>(Istream& is, LList<LListBase, T>& L)",
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        is.fatalCheck(" operator>>(Istream& is, LList<LListBase, T>& L)");

        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            T element;
            is >> element;
            L.append(element);

            is >> lastToken;
            is.fatalCheck(" operator>>(Istream& is, LList<LListBase, T>& L)");
        }
    }
    else
    {
        FatalIOErrorIn(" operator>>(Istream& is, LList<LListBase, T>& L)", is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    // Check state of IOstream
    is.fatalCheck(" operator>>(Istream& is, LList<LListBase, T>& L)");

    return is;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Ostream& operator<<(Ostream& os, const LList<LListBase, T>& ll)
{
    // Write size of LList
    os << nl << ll.size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST << nl;

    // Write LList contents
    for
    (
        typename LList<LListBase, T>::const_iterator iter = ll.begin();
        iter != ll.end();
        ++iter
    )
    {
        os << iter() << nl;
    }

    // Write end of contents
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const LList&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
