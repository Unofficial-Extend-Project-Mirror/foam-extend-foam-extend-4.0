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

#include "directTopoChange.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Foam::directTopoChange::reorder
(
    const labelList& oldToNew,
    DynamicList<T>& lst
)
{
    // Create copy
    DynamicList<T> oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI] = oldLst[elemI];
        }
    }
}


template <class T>
void Foam::directTopoChange::reorder
(
    const labelList& oldToNew,
    List<DynamicList<T> >& lst
)
{
    // Create copy
    List<DynamicList<T> > oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI].transfer(oldLst[elemI]);
        }
    }
}


template <class T>
void Foam::directTopoChange::renumberKey
(
    const labelList& oldToNew,
    Map<T>& elems
)
{
    Map<T> newElems(elems.size());

    forAllConstIter(typename Map<T>, elems, iter)
    {
        label newElem = oldToNew[iter.key()];

        if (newElem >= 0)
        {
            newElems.insert(newElem, iter());
        }
    }

    elems.transfer(newElems);
}


// ************************************************************************* //
