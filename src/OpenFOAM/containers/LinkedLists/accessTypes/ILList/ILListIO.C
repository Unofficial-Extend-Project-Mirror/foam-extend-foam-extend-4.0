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

#include "ILList.H"
#include "Istream.H"
#include "INew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
void ILList<LListBase, T>::read(Istream& is, const INew& inewt)
{
    is.fatalCheck("operator>>(Istream& is, ILList<LListBase, T>& L)");

    token firstToken(is);

    is.fatalCheck
    (
        "operator>>(Istream& is, ILList<LListBase, T>& L) : reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char listDelimiter = is.readBeginList("ILList<LListBase, T>");

        if (s)
        {
            if (listDelimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; i++)
                {
                    append(inewt(is).ptr());
                
                    is.fatalCheck
                    (
                        "operator>>(Istream& is, ILList<LListBase, T>& L) : "
                        "reading entry"
                    );
                }
            }
            else
            {
                T* tPtr = inewt(is).ptr();
                append(tPtr);

                is.fatalCheck
                (
                    "operator>>(Istream& is, ILList<LListBase, T>& L) : "
                    "reading entry"
                );
                
                for (label i=1; i<s; i++)
                {
                    append(new T(*tPtr));
                }
            }
        }

        // Read end of contents
        is.readEndList("ILList<LListBase, T>");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "operator>>(Istream& is, ILList<LListBase, T>& L)",
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        is.fatalCheck("operator>>(Istream& is, ILList<LListBase, T>& L)");

        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            append(inewt(is).ptr());

            is >> lastToken;
            is.fatalCheck("operator>>(Istream& is, ILList<LListBase, T>& L)");
        }
    }
    else
    {
        FatalIOErrorIn("operator>>(Istream& is, ILList<LListBase, T>& L)", is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck("operator>>(Istream& is, ILList<LListBase, T>& L)");
}


//- Construct from Istream using given Istream constructor class
template<class LListBase, class T>
template<class INew>
ILList<LListBase, T>::ILList(Istream& is, const INew& inewt)
{
    read(is, inewt);
}


// Construct from Istream
template<class LListBase, class T>
ILList<LListBase, T>::ILList(Istream& is)
{
    read(is, INew<T>());
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Istream& operator>>(Istream& is, ILList<LListBase, T>& L)
{
    L.clear();
    L.read(is, INew<T>());

    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
