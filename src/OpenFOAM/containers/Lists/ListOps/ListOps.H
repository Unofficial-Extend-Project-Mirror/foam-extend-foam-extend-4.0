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

InNamspace
    Foam

Description
    Various functions to operate on Lists.

SourceFiles
    ListOps.C
    ListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ListOps_H
#define ListOps_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Renumber the values (not the indices) of a list. List elements <= 0 are
//  left as is.
template<class List>
List renumber(const labelList& oldToNew, const List&);

//- Inplace renumber the values of a list. List elements <= 0 are
//  left as is.
template<class List>
void inplaceRenumber(const labelList& oldToNew, List&);


//- Reorder the elements (indices, not values) of a list.
//  List elements <= 0 are left as is.
template<class List>
List reorder(const labelList& oldToNew, const List&);

//- Inplace reorder the elements of a list.
//  List elements <= 0 are left as is.
template<class List>
void inplaceReorder(const labelList& oldToNew, List&);


// Variants to work with iterators and sparse tables. Need to have iterators
// and insert()

//- Map values. Do not map negative values.
template<class Container>
void inplaceMapValue(const labelList& oldToNew, Container&);
//- Recreate with mapped keys. Remove elements with negative key.
template<class Container>
void inplaceMapKey(const labelList& oldToNew, Container&);


//- Extract elements of List whose region is certain value. Use e.g.
//  to extract all selected elements:
//  subset<boolList, labelList>(selectedElems, true, lst);
template<class T, class List>
List subset(const UList<T>& regions, const T& region, const List&);

//- Inplace extract elements of List whose region is certain value. Use e.g.
//  to extract all selected elements:
//  inplaceSubset<boolList, labelList>(selectedElems, true, lst);
template<class T, class List>
void inplaceSubset(const UList<T>& regions, const T& region, List&);

//- Invert one-to-one map. Unmapped elements will be -1.
labelList invert(const label len, const labelList& oldToNew);

//- Invert one-to-many map. Unmapped elements will be size 0.
labelListList invertOneToMany(const label len, const labelList&);

//- Invert many-to-many. Input and output types need to be inherited
//  from List. E.g. faces to pointFaces.
template<class InList, class OutList>
void invertManyToMany(const label len, const UList<InList>&, List<OutList>&);

template<class InList, class OutList>
List<OutList> invertManyToMany(const label len, const UList<InList>& in)
{
    List<OutList> out;
    invertManyToMany<InList,OutList>(len, in, out);
    return out;
}

//- Create identity map (map[i] == i) of given length
labelList identity(const label len);

//- Find first occurence of given element and return index,
//  return -1 if not found. Linear search.
template<class List>
label findIndex
(
    const List&,
    typename List::const_reference,
    const label start=0
);

//- Find all occurences of given element. Linear search.
template<class List>
labelList findIndices
(
    const List&,
    typename List::const_reference,
    const label start=0
);

//- Opposite of findIndices: set values at indices to given value
template<class List>
void setValues
(
    List&,
    const labelList& indices,
    typename List::const_reference
);

//- Opposite of findIndices: set values at indices to given value
template<class List>
List createWithValues
(
    const label sz,
    typename List::const_reference initValue,
    const labelList& indices,
    typename List::const_reference setValue
);

//- Find index of max element (and larger than given element).
//  return -1 if not found. Linear search.
template<class List>
label findMax(const List&, const label start=0);


//- Find index of min element (and less than given element).
//  return -1 if not found. Linear search.
template<class List>
label findMin(const List&, const label start=0);


//- Find first occurence of given element in sorted list and return index,
//  return -1 if not found. Binary search.
template<class List>
label findSortedIndex
(
    const List&,
    typename List::const_reference,
    const label start=0
);


//- Find last element < given value in sorted list and return index,
//  return -1 if not found. Binary search.
template<class List>
label findLower
(
    const List&,
    typename List::const_reference,
    const label start=0
);


//- To construct a List from a C array. Has extra Container type
//  to initialise e.g. wordList from arrays of char*.
template<class Container, class T, int nRows>
List<Container> initList(const T[nRows]);


//- To construct a (square) ListList from a C array. Has extra Container type
//  to initialise e.g. faceList from arrays of labels.
template<class Container, class T, int nRows, int nColumns>
List<Container> initListList(const T[nRows][nColumns]);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

