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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList()
:
    List<T>(SizeInc),
    nextFree_(0)
{}


//- Construct given size
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const label s
)
:
    List<T>(s),
    nextFree_(0)
{}


//- Construct given size
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const UList<T>& s
)
:
    List<T>(s),
    nextFree_(s.size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::allocSize()
const
{
    return List<T>::size();
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::size()
const
{
    return nextFree_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s
)
{
    if (s < nextFree_)
    {
        nextFree_ = s;
    }
    else
    {
        List<T>::setSize(s);
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s,
    const T& t
)
{
    if (s < nextFree_)
    {
        nextFree_ = s;
    }
    else
    {
        List<T>::setSize(s, t);
        nextFree_ = s;
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clear()
{
    nextFree_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::shrink()
{
    List<T>::setSize(nextFree_);
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer(List<T>& l)
{
    if (l.size() < nextFree_)
    {
        FatalErrorIn
        (
            "void DynamicList<T, SizeInc, SizeMult"
            ", SizeDiv>::transfer(List<T>&)"
        )   << "Cannot replace the underlying storage of this DynamicList"
            << " of which " << nextFree_ << " elements are used" << nl
            << "with a List of size " << l.size() << abort(FatalError);
    }
    else
    {
        List<T>::transfer(l);   // take over storage
        l.clear();              // set nextFree of l to 0
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer
(
    DynamicList<T, SizeInc, SizeMult, SizeDiv>& l
)
{
    nextFree_ = l.nextFree_;    // take over used size
    // this should disable shrink() in List<T>::transfer - but somehow doesn't
    l.nextFree_ = l.allocSize();
    List<T>::transfer(l);       // take over storage
    l.clear();                  // set nextFree of l to 0
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const T& e
)
{
    nextFree_++;

    if (nextFree_ > List<T>::size())
    {
        List<T>::setSize
        (
            max
            (
                nextFree_,
                label(SizeMult*List<T>::size()/SizeDiv + SizeInc)
            )
        );
    }

    this->operator[](nextFree_ - 1) = e;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()
{
    if (nextFree_ == 0)
    {
        FatalErrorIn
        (
            "void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()"
        )   << "List is empty" << abort(FatalError);
    }

    return List<T>::operator[](--nextFree_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T& Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator()
(
    const label i
)
{
    nextFree_ = max(nextFree_, i + 1);

    if (nextFree_ > List<T>::size())
    {
        List<T>::setSize
        (
            max
            (
                nextFree_,
                label(SizeMult*List<T>::size()/SizeDiv + SizeInc)
            )
        );
    }

    return this->operator[](i);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const T& t
)
{
    List<T>::operator=(t);
    nextFree_ = List<T>::size();
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const List<T>& l
)
{
    List<T>::operator=(l);
    nextFree_ = l.size();
}


// ************************************************************************* //
