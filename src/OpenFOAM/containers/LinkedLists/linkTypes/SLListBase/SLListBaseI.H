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
    Base singly-linked list.

\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

inline SLListBase::link::link()
:
    next_(0)
{}


inline SLListBase::link::link(link* p)
:
    next_(p)
{}


inline SLListBase::SLListBase()
:
    last_(0),
    nElmts_(0)
{}


inline SLListBase::SLListBase(link* a)
:
    last_(a->next_ = a),
    nElmts_(1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline SLListBase::~SLListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label SLListBase::size() const
{
    return nElmts_;
}


inline SLListBase::link* SLListBase::first()
{
    if (!nElmts_)
    {
        FatalErrorIn("SLListBase::first()")
            << "list is empty"
            << abort(FatalError);
    }
    return last_->next_;
}


inline const SLListBase::link* SLListBase::first() const
{
    if (!nElmts_)
    {
        FatalErrorIn("SLListBase::first() const")
            << "list is empty"
            << abort(FatalError);
    }
    return last_->next_;
}


inline SLListBase::link* SLListBase::last()
{
    if (!nElmts_)
    {
        FatalErrorIn("SLListBase::last()")
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline const SLListBase::link* SLListBase::last() const
{
    if (!nElmts_)
    {
        FatalErrorIn("SLListBase::last() const")
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline void SLListBase::clear()
{
    nElmts_ = 0;
    last_ = 0;
}


inline SLListBase::link* SLListBase::remove(SLListBase::iterator& it)
{
    return remove(it.curElmt_);
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

inline SLListBase::iterator::iterator(SLListBase& s, link* elmt)
:
    curList_(s),
    curElmt_(elmt),
    curLink_(*curElmt_)
{}


inline SLListBase::iterator::iterator(SLListBase& s)
:
    curList_(s),
    curElmt_(NULL),
    curLink_()
{}


inline void SLListBase::iterator::operator=(const iterator& iter)
{
    curElmt_ = iter.curElmt_;
    curLink_ = iter.curLink_;
}


inline bool SLListBase::iterator::operator==(const iterator& iter) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool SLListBase::iterator::operator!=(const iterator& iter) const
{
    return curElmt_ != iter.curElmt_;
}


inline SLListBase::link& SLListBase::iterator::operator*()
{
    return *curElmt_;
}


inline SLListBase::iterator& SLListBase::iterator::operator++()
{
    if (curElmt_ == curList_.last_ || curList_.last_ == 0)
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


inline SLListBase::iterator SLListBase::iterator::operator++(int)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


inline SLListBase::iterator SLListBase::begin()
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


inline const SLListBase::iterator& SLListBase::end()
{
    return endIter;
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

inline SLListBase::const_iterator::const_iterator
(
    const SLListBase& s,
    const link* elmt
)
:
    curList_(s),
    curElmt_(elmt)
{}


inline SLListBase::const_iterator::const_iterator(const iterator& iter)
:
    curList_(iter.curList_),
    curElmt_(iter.curElmt_)
{}


inline void SLListBase::const_iterator::operator=(const const_iterator& iter)
{
    curElmt_ = iter.curElmt_;
}


inline bool SLListBase::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool SLListBase::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return curElmt_ != iter.curElmt_;
}


inline const SLListBase::link& SLListBase::const_iterator::operator*()
{
    return *curElmt_;
}


inline SLListBase::const_iterator& SLListBase::const_iterator::operator++()
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


inline SLListBase::const_iterator SLListBase::const_iterator::operator++(int)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


inline SLListBase::const_iterator SLListBase::begin() const
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


inline const SLListBase::const_iterator& SLListBase::end() const
{
    return endConstIter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
