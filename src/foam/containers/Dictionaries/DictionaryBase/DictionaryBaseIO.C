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

Description
    Reads the data description and data portions of a DictionaryBase File.

\*---------------------------------------------------------------------------*/

#include "DictionaryBase.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

template<class IDLListType, class T>
Ostream& operator<<(Ostream& os, const DictionaryBase<IDLListType, T>& dict)
{
    for
    (
        typename IDLListType::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        os << *iter;

        // Check stream before going to next entry.
        if (!os.good())
        {
            Info
                << "operator<<(Ostream& os, const DictionaryBase&) : "
                << "Can't write entry for DictionaryBase"
                << endl;

            return os;
        }
    }

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
