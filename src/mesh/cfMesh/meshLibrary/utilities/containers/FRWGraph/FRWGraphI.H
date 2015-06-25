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

\*---------------------------------------------------------------------------*/

template<class T, Foam::label width>
void Foam::FRWGraph<T, width>::checkIndex(const label i, const label j) const
{
    if( (i < 0) || (i >= nRows_) )
    {
        FatalErrorIn
        (
            "void Foam::FRWGraph<T,width>::"
            "checkIndex(const label i, const label j) const"
        ) << "Row index " << i
            << " is not in range " << 0
            << " and " << nRows_ << abort(FatalError);
    }
    
    if( (j < 0) || (j >= width) )
        FatalErrorIn
        (
            "void Foam::FRWGraph<T,width>::"
            "checkIndex(const label i, const label j) const"
        ) << "Column index " << j
            << " is not in range " << 0
            << " and " << width << abort(FatalError);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T, Foam::label width>
inline Foam::FRWGraph<T,width>::FRWGraph()
:
    data_(),
    nRows_(0)
{}

//- Construct given size
template<class T, Foam::label width>
inline Foam::FRWGraph<T,width>::FRWGraph
(
    const label s
)
:
    data_(s * width),
    nRows_(s)
{}


//- Construct given size
template<class T, Foam::label width>
inline Foam::FRWGraph<T,width>::FRWGraph
(
    const label s,
    const T& t
)
:
    data_(s * width, t),
    nRows_(s)
{}

template<class T, Foam::label width>
inline Foam::FRWGraph<T,width>::FRWGraph
(
    const FRWGraph<T,width>& ol
)
:
    data_(ol.data_),
    nRows_(ol.nRows_)
{}

template<class T, Foam::label width>
inline Foam::FRWGraph<T,width>::~FRWGraph()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, Foam::label width>
inline Foam::label Foam::FRWGraph<T,width>::size() const
{
    return nRows_;
}

template<class T, Foam::label width>
inline Foam::label Foam::FRWGraph<T,width>::sizeOfRow(const label rowI) const
{
    return width;
}

template<class T, Foam::label width>
inline void Foam::FRWGraph<T,width>::setSize(const label i)
{
    data_.setSize(i * width);
    nRows_ = i;
}

template<class T, Foam::label width>
inline void Foam::FRWGraph<T,width>::clear()
{
    data_.clear();
    nRows_ = 0;
}

template<class T, Foam::label width>
inline void Foam::FRWGraph<T,width>::appendFixedList
(
    const FixedList<T, width>& l
)
{
    forAll(l, elI)
        data_.append(l[elI]);
    ++nRows_;
}

template<class T, Foam::label width>
inline void Foam::FRWGraph<T,width>::setRow
(
    const label rowI,
    const FixedList<T, width>& l
)
{
    const label start = rowI * width;
    forAll(l, elI)
        data_[start+elI] = l[elI];
}

template<class T, Foam::label width>
inline bool Foam::FRWGraph<T,width>::contains
(
    const label rowI,
    const T& e
) const
{
    const label start = rowI * width;
    for(register label i=0;i<width;++i)
        if( data_[start+i] == e )
            return true;
    
    return false;
}

template<class T, Foam::label width>
inline Foam::label Foam::FRWGraph<T,width>::containsAtPosition
(
    const label rowI,
    const T& e
) const
{
    const label start = rowI * width;
    for(register label i=0;i<width;++i)
        if( data_[start+i] == e )
            return i;
    
    return -1;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
template<class T, Foam::label width>
inline const T& Foam::FRWGraph<T,width>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif
        
    return data_[i * width + j];
}

template<class T, Foam::label width>
inline T& Foam::FRWGraph<T,width>::operator()
(
    const label i, const label j
)
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif
        
    return data_[i * width + j];
}

template<class T, Foam::label width>
inline void Foam::FRWGraph<T,width>::operator=
(
    const FRWGraph<T, width>& l
)
{
    data_ = l.data_;
    nRows_ = l.nRows_;
}


// ************************************************************************* //
