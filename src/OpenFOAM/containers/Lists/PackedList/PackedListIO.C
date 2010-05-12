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

#include "PackedList.H"
#include "Ostream.H"
#include "token.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

template<class T>
void PackedList<T>::writeEntry(Ostream& os) const
{
    if
    (
        size()
     && token::compound::isCompound
        (
            "List<" + word(pTraits<T>::typeName) + '>'
        )
    )
    {
        os  << word("List<" + word(pTraits<T>::typeName) + '>') << " ";
    }
    
    os << *this;
}


template<class T>
void PackedList<T>::writeEntry(const word& keyword, Ostream& os) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T>
Ostream& operator<<(Ostream& os, const PackedList<T>& L)
{
    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        bool uniform = false;

        if (L.size() > 1 && contiguous<T>())
        {
            uniform = true;

            unsigned int L0 = L.get(0);

            forAll(L, i)
            {
                if (L.get(i) != L0)
                {
                    uniform = false;
                    break;
                }
            }
        }

        if (uniform)
        {
            // Write size of list and start contents delimiter
            os << L.size() << token::BEGIN_BLOCK;

            // Write list contents
            os << L.get(0);

            // Write end of contents delimiter
            os << token::END_BLOCK;
        }
        else if (L.size() < 11 && contiguous<T>())
        {
            // Write size of list and start contents delimiter
            os << L.size() << token::BEGIN_LIST;

            // Write list contents
            forAll(L, i)
            {
                if (i > 0) os << token::SPACE;
                os << L.get(i);
            }

            // Write end of contents delimiter
            os << token::END_LIST;
        }
        else
        {
            // Write size of list and start contents delimiter
            os << nl << L.size() << nl << token::BEGIN_LIST;

            // Write list contents
            forAll(L, i)
            {
                os << nl << L.get(i);
            }

            // Write end of contents delimiter
            os << nl << token::END_LIST << nl;
        }
    }
    else
    {
        os << nl << L.size() << nl;
        if (L.size())
        {
            os.write((const char*)L.v_, L.size()*sizeof(T));
        }
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const PackedList&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
