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

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

inline DLListBase::link::link()
:
    prev_(0),
    next_(0)
{}


inline DLListBase::DLListBase()
:
    first_(0),
    last_(0),
    nElmts_(0)
{}


inline DLListBase::DLListBase(link* a)
:
    first_(a),
    last_(a),
    nElmts_(1)
{
    a->prev_ = a;
    a->next_ = a;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline DLListBase::~DLListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool DLListBase::link::registered() const
{
    return prev_ != 0 && next_ != 0;
}


inline void DLListBase::link::deregister()
{
    prev_ = 0;
    next_ = 0;
}


inline label DLListBase::size() const
{
    return nElmts_;
}


inline DLListBase::link* DLListBase::first()
{
    if (!nElmts_)
    {
        FatalErrorIn("DLListBase::first()")
            << "list is empty"
            << abort(FatalError);
    }
    return first_;
}


inline const DLListBase::link* DLListBase::first() const
{
    if (!nElmts_)
    {
        FatalErrorIn("DLListBase::first() const")
            << "list is empty"
            << abort(FatalError);
    }
    return first_;
}


inline DLListBase::link* DLListBase::last()
{
    if (!nElmts_)
    {
        FatalErrorIn("DLListBase::last()")
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline const DLListBase::link* DLListBase::last() const
{
    if (!nElmts_)
    {
        FatalErrorIn("DLListBase::last() const")
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline void DLListBase::clear()
{
    nElmts_ = 0;
    first_ = 0;
    last_ = 0;
}


inline DLListBase::link* DLListBase::remove(DLListBase::iterator& it)
{
    return remove(it.curElmt_);
}


inline DLListBase::link* DLListBase::replace
(
    DLListBase::iterator& oldIter,
    DLListBase::link* newLink
)
{
    return replace(oldIter.curElmt_, newLink);
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

inline DLListBase::iterator::iterator(DLListBase& s, link* elmt)
:
    curList_(s),
    curElmt_(elmt),
    curLink_(*curElmt_)
{}


inline DLListBase::iterator::iterator(DLListBase& s)
:
    curList_(s),
    curElmt_(NULL),
    curLink_()
{}


inline void DLListBase::iterator::operator=(const iterator& iter)
{
    curElmt_ = iter.curElmt_;
    curLink_ = iter.curLink_;
}


inline bool DLListBase::iterator::operator==(const iterator& iter) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool DLListBase::iterator::operator!=(const iterator& iter) const
{
    return curElmt_ != iter.curElmt_;
}


inline DLListBase::link& DLListBase::iterator::operator*()
{
    return *curElmt_;
}


inline DLListBase::iterator& DLListBase::iterator::operator++()
{
    // Check if the curElmt_ is the last element (if it points to itself)
    // or if the list is empty because the last element may have been removed
    if (curLink_.next_ == curElmt_ || curList_.last_ == 0)
    {
        curElmt_ = 0;
    }
    else
    {
        curElmt_ = curLink_.next_;
        curLink_ = *curElmt_;
    }

    return *this;
}


inline DLListBase::iterator DLListBase::iterator::operator++(int)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


inline DLListBase::iterator DLListBase::begin()
{
    if (size())
    {
        return iterator(*this, first());
    }
    else
    {
        return endIter;
    }
}


inline const DLListBase::iterator& DLListBase::end()
{
    return endIter;
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

inline DLListBase::const_iterator::const_iterator
(
    const DLListBase& s,
    const link* elmt
)
:
    curList_(s),
    curElmt_(elmt)
{}


inline DLListBase::const_iterator::const_iterator(const iterator& iter)
:
    curList_(iter.curList_),
    curElmt_(iter.curElmt_)
{}


inline void DLListBase::const_iterator::operator=(const const_iterator& iter)
{
    curElmt_ = iter.curElmt_;
}


inline bool DLListBase::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool DLListBase::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return curElmt_ != iter.curElmt_;
}


inline const DLListBase::link& DLListBase::const_iterator::operator*()
{
    return *curElmt_;
}


inline DLListBase::const_iterator& DLListBase::const_iterator::operator++()
{
    if (curElmt_ == curList_.last_)
    {
        curElmt_ = 0;
    }
    else
    {
        curElmt_ = curElmt_->next_;
    }

    return *this;
}


inline DLListBase::const_iterator DLListBase::const_iterator::operator++(int)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


inline DLListBase::const_iterator DLListBase::begin() const
{
    if (size())
    {
        return const_iterator(*this, first());
    }
    else
    {
        return endConstIter;
    }
}


inline const DLListBase::const_iterator& DLListBase::end() const
{
    return endConstIter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
