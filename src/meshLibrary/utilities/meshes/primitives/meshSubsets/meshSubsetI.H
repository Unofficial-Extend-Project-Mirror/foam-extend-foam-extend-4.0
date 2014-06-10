/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "meshSubset.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline meshSubset::meshSubset()
:
    name_(),
    type_(UNKNOWN),
    data_()
{}

inline meshSubset::meshSubset
(
    const word& name,
    const meshSubset::subsetType_& t
)
:
    name_(name),
    type_(t),
    data_()
{}

template<class ListType>
inline meshSubset::meshSubset
(
    const word& name,
    const meshSubset::subsetType_& type,
    const ListType& elements
)
:
    name_(name),
    type_(type),
    data_()
{
    forAll(elements, i)
        data_.insert(elements[i]);
}

inline meshSubset::meshSubset(const meshSubset& ms)
:
    name_(ms.name_),
    type_(ms.type_),
    data_()
{
    forAllConstIter(HashSet<label>, ms.data_, it)
        data_.insert(it.key());
}

inline meshSubset::meshSubset(Istream& is)
:
    name_(),
    type_(UNKNOWN),
    data_()
{
    is >> *this;
}

inline meshSubset::~meshSubset()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& meshSubset::name() const
{
    return name_;
}

inline label meshSubset::type() const
{
    return type_;
}

template<class ListType>
inline void meshSubset::containedElements(ListType& l) const
{
    l.setSize(data_.size());

    label counter(0);
    forAllConstIter(HashSet<label>, data_, it)
        l[counter++] = it.key();
}

inline void meshSubset::addElement(const label elmt)
{
    data_.insert(elmt);
}

inline void meshSubset::removeElement(const label elmt)
{
    data_.erase(elmt);
}

template<class ListType>
inline void meshSubset::updateSubset(const ListType& newLabels)
{
    HashSet<label> newData;

    forAllConstIter(HashSet<label>, data_, it)
    {
        if( newLabels[it.key()] < 0 )
            continue;

        newData.insert(newLabels[it.key()]);
    }

    data_.clear();
    forAllConstIter(HashSet<label>, newData, it)
        data_.insert(it.key());
}

inline void meshSubset::updateSubset(const VRWGraph& newLabels)
{
    HashSet<label> newData;

    forAllConstIter(HashSet<label>, data_, it)
    {
        forAllRow(newLabels, it.key(), i)
            newData.insert(newLabels(it.key(), i));
    }

    data_.clear();
    forAllConstIter(HashSet<label>, newData, it)
        data_.insert(it.key());
}

inline bool meshSubset::contains(const label elmt) const
{
    return data_.found(elmt);
}

inline void meshSubset::operator=(const meshSubset& ms)
{
    name_ = ms.name_;
    type_ = ms.type_;

    data_.clear();
    forAllConstIter(HashSet<label>, ms.data_, it)
        data_.insert(it.key());
}

inline bool meshSubset::operator==(const meshSubset& ms) const
{
    if( ms.name_ != name_ )
        return false;
    if( ms.type_ != ms.type_ )
        return false;

    forAllConstIter(HashSet<label>, data_, it)
        if( !ms.data_.found(it.key()) )
            return false;

    return true;
}

inline bool meshSubset::operator!=(const meshSubset& ms) const
{
    return !operator==(ms);
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const meshSubset& sel)
{
    os.check("inline Ostream& operator<<(Ostream&, const meshSubset&)");

    os << sel.name_ << nl << sel.type_ << nl << sel.data_;

    return os;
}

inline Istream& operator>>(Istream& is, meshSubset& sel)
{
    is.check("friend Istream& operator>>(Istream&, meshSubset&)");

    is >> sel.name_ >> sel.type_ >> sel.data_;

    is.check("friend Istream& operator>>(Istream&, meshSubset&)");

    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
