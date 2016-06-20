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

Class
    Foam::DLPtrList

Description
    Non-intrusive doubly-linked pointer list.

\*---------------------------------------------------------------------------*/

#ifndef DLPtrList_H
#define DLPtrList_H

#include "LPtrList.H"
#include "DLListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DLPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class DLPtrList
:
    public LPtrList<DLListBase, T>
{

public:

    // Constructors

        //- Null construct
        DLPtrList()
        {}

        //- Construct given initial T
        DLPtrList(T a)
        :
            LPtrList<DLListBase, T>(a)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        DLPtrList(Istream& is, const INew& inewt)
        :
            LPtrList<DLListBase, T>(is, inewt)
        {}

        //- Construct from Istream
        DLPtrList(Istream& is)
        :
            LPtrList<DLListBase, T>(is)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
