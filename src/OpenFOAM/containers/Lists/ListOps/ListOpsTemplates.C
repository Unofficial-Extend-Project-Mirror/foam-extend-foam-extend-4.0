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

#include "ListOps.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class List>
List Foam::renumber
(
    const labelList& oldToNew,
    const List& lst
)
{
    // Create copy
    List newLst(lst.size());

    forAll(lst, elemI)
    {
        if (lst[elemI] >= 0)
        {
            newLst[elemI] = oldToNew[lst[elemI]];
        }
    }

    return newLst;
}


template<class List>
void Foam::inplaceRenumber
(
    const labelList& oldToNew,
    List& lst
)
{
    forAll(lst, elemI)
    {
        if (lst[elemI] >= 0)
        {
            lst[elemI] = oldToNew[lst[elemI]];
        }
    }
}


template<class List>
List Foam::reorder
(
    const labelList& oldToNew,
    const List& lst
)
{
    // Create copy
    List newLst(lst.size());

    forAll(lst, elemI)
    {
        if (oldToNew[elemI] >= 0)
        {
            newLst[oldToNew[elemI]] = lst[elemI];
        }
        else
        {
            newLst[elemI] = lst[elemI];
        }
    }
    return newLst;
}


template<class List>
void Foam::inplaceReorder
(
    const labelList& oldToNew,
    List& lst
)
{
    // Create copy
    List newLst(lst.size());

    forAll(lst, elemI)
    {
        if (oldToNew[elemI] >= 0)
        {
            newLst[oldToNew[elemI]] = lst[elemI];
        }
        else
        {
            newLst[elemI] = lst[elemI];
        }
    }

    lst.transfer(newLst);
}


template<class Container>
void Foam::inplaceMapValue
(
    const labelList& oldToNew,
    Container& lst
)
{
    for
    (
        typename Container::iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        if (iter() >= 0)
        {
            iter() = oldToNew[iter()];
        }
    }
}


template<class Container>
void Foam::inplaceMapKey
(
    const labelList& oldToNew,
    Container& lst
)
{
    Container newLst(lst);

    for
    (
        typename Container::iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        if (iter.key() >= 0)
        {
            newLst.insert(oldToNew[iter.key()], iter());
        }
    }
    
    lst.transfer(newLst);
}


template<class T, class List>
List Foam::subset(const UList<T>& regions, const T& region, const List& lst)
{
    if (regions.size() < lst.size())
    {
        FatalErrorIn("subset(const UList<T>&, const T&, const List&)")
            << "Regions is of size " << regions.size()
            << "; list it is supposed to index is of size " << lst.size()
            << abort(FatalError);
    }

    List newLst(lst.size());

    label nElem = 0;
    forAll(lst, elemI)
    {
        if (regions[elemI] == region)
        {
            newLst[nElem++] = lst[elemI];
        }
    }
    newLst.setSize(nElem);

    return newLst;
}


template<class T, class List>
void Foam::inplaceSubset(const UList<T>& regions, const T& region, List& lst)
{
    if (regions.size() < lst.size())
    {
        FatalErrorIn("inplaceSubset(const UList<T>&, const T&, List&)")
            << "Regions is of size " << regions.size()
            << "; list it is supposed to index is of size " << lst.size()
            << abort(FatalError);
    }

    label nElem = 0;
    forAll(lst, elemI)
    {
        if (regions[elemI] == region)
        {
            if (nElem != elemI)
            {
                lst[nElem] = lst[elemI];
            }
            ++nElem;
        }
    }

    lst.setSize(nElem);
}


// As clarification coded as inversion from pointEdges to edges but completely
// general.
template<class InList, class OutList>
void Foam::invertManyToMany
(
    const label nEdges,
    const UList<InList>& pointEdges,
    List<OutList>& edges
)
{
    // Number of points per edge
    labelList nPointsPerEdge(nEdges, 0);

    forAll(pointEdges, pointI)
    {
        const InList& pEdges = pointEdges[pointI];

        forAll(pEdges, j)
        {
            nPointsPerEdge[pEdges[j]]++;
        }
    }

    // Size edges
    edges.setSize(nEdges);

    forAll(nPointsPerEdge, edgeI)
    {
        edges[edgeI].setSize(nPointsPerEdge[edgeI]);
    }
    nPointsPerEdge = 0;

    // Fill edges
    forAll(pointEdges, pointI)
    {
        const InList& pEdges = pointEdges[pointI];

        forAll(pEdges, j)
        {
            label edgeI = pEdges[j];

            edges[edgeI][nPointsPerEdge[edgeI]++] = pointI;
        }
    }
}


template<class List>
Foam::label Foam::findIndex
(
    const List& l,
    typename List::const_reference t,
    const label start
)
{
    label index = -1;

    for (label i = start; i < l.size(); i++)
    {
        if (l[i] == t)
        {
            index = i;
            break;
        }
    }

    return index;
}


template<class List>
Foam::labelList Foam::findIndices
(
    const List& l,
    typename List::const_reference t,
    const label start
)
{
    // Count occurrences
    label n = 0;

    for (label i = start; i < l.size(); i++)
    {
        if (l[i] == t)
        {
            n++;
        }
    }

    // Create and fill
    labelList indices(n);
    n = 0;

    for (label i = start; i < l.size(); i++)
    {
        if (l[i] == t)
        {
            indices[n++] = i;
        }
    }

    return indices;
}


template<class List>
void Foam::setValues
(
    List& l,
    const labelList& indices,
    typename List::const_reference t
)
{
    forAll(indices, i)
    {
        l[indices[i]] = t;
    }
}


template<class List>
List Foam::createWithValues
(
    const label sz,
    const typename List::const_reference initValue,
    const labelList& indices,
    typename List::const_reference setValue
)
{
    List l(sz, initValue);
    setValues(l, indices, setValue);
    return l;
}


template<class List>
Foam::label Foam::findMax(const List& l, const label start)
{
    if (start >= l.size())
    {
        return -1;
    }

    label index = start;

    for (label i = start+1; i < l.size(); i++)
    {
        if (l[i] > l[index])
        {
            index = i;
        }
    }

    return index;
}


template<class List>
Foam::label Foam::findMin(const List& l, const label start)
{
    if (start >= l.size())
    {
        return -1;
    }

    label index = start;

    for (label i = start+1; i < l.size(); i++)
    {
        if (l[i] < l[index])
        {
            index = i;
        }
    }

    return index;
}


template<class List>
Foam::label Foam::findSortedIndex
(
    const List& l,
    typename List::const_reference t,
    const label start
)
{
    if (start >= l.size())
    {
        return -1;
    }

    label low = start;
    label high = l.size() - 1;

    while (low <= high)
    {
        label mid = (low + high)/2;

        if (t < l[mid])
        {
            high = mid - 1;
        }
        else if (t > l[mid])
        {
            low = mid + 1;
        }
        else
        {
            return mid;
        }
    }

    return -1;
}


template<class List>
Foam::label Foam::findLower
(
    const List& l,
    typename List::const_reference t,
    const label start
)
{
    if (start >= l.size())
    {
        return -1;
    }

    label low = start;
    label high = l.size() - 1;

    while ((high - low) > 1)
    {
        label mid = (low + high)/2;

        if (l[mid] < t)
        {
            low = mid;
        }
        else
        {
            high = mid;
        }
    }

    if (l[high] < t)
    {
        return high;
    }
    else
    {
        if (l[low] < t)
        {
            return low;
        }
        else
        {
            return -1;
        }
    }
}


template<class Container, class T, int nRows>
Foam::List<Container> Foam::initList(const T elems[nRows])
{
    List<Container> faces(nRows);

    forAll(faces, faceI)
    {
        faces[faceI] = Container(elems[faceI]);
    }
    return faces;
}


template<class Container, class T, int nRows, int nColumns>
Foam::List<Container> Foam::initListList(const T elems[nRows][nColumns])
{
    List<Container> faces(nRows);

    Container f(nColumns);
    forAll(faces, faceI)
    {
        forAll(f, i)
        {
            f[i] = elems[faceI][i];
        }
        faces[faceI] = f;
    }
    return faces;
}


// ************************************************************************* //
