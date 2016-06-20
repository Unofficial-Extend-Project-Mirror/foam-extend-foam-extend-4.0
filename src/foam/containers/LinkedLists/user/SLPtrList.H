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
    Foam::SLPtrList

Description
    Non-intrusive singly-linked pointer list.

\*---------------------------------------------------------------------------*/

#ifndef SLPtrList_H
#define SLPtrList_H

#include "LPtrList.H"
#include "SLListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SLPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SLPtrList
:
    public LPtrList<SLListBase, T>
{

public:

    // Constructors

        //- Null construct
        SLPtrList()
        {}

        //- Construct given initial T*
        SLPtrList(T* a)
        :
            LPtrList<SLListBase, T>(a)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        SLPtrList(Istream& is, const INew& inewt)
        :
            LPtrList<SLListBase, T>(is, inewt)
        {}

        //- Construct from Istream
        SLPtrList(Istream& is)
        :
            LPtrList<SLListBase, T>(is)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
