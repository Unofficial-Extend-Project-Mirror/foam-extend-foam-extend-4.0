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

#include "contiguous.H"

template<class T, Foam::label Offset>
void Foam::LongList<T, Offset>::checkIndex(const label i) const
{
    if( (i < 0) || (i >= nextFree_) )
    {
        FatalErrorIn
        (
            "void Foam::LongList<T, label, Offset>::"
            "checkIndex(const label i) const"
        ) << "Index " << Foam::label(i) << " is not in range " << Foam::label(0)
            << " and " << Foam::label(nextFree_) << abort(FatalError);
    }
}

template<class T, Foam::label Offset>
void Foam::LongList<T, Offset>::initializeParameters()
{
    unsigned int t = sizeof(T);
    label it(0);
    
    while( t > 1 )
    {
        t >>= 1;
        ++it;
    }
    
    shift_ = Foam::max(10, Offset - it);
    mask_ = 1<<shift_;
    mask_ -= 1;
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::allocateSize(const label s)
{
    if( s == 0 )
    {
        clearOut();
        return;
    }
    
    const label numblock1 = ((s-1)>>shift_) + 1;
    const label blockSize = 1<<shift_;
   
    if( numblock1 < numBlocks_ )
    {
        for(register label i=numblock1;i<numBlocks_;++i)
            delete [] dataPtr_[i];
    }
    else if( numblock1 > numBlocks_ )
    {
        if( numblock1 >= numAllocatedBlocks_ )
        {
            do
            {
                numAllocatedBlocks_ += 64;
            } while( numblock1 > numAllocatedBlocks_ );
            
            T** dataptr1 = new T*[numAllocatedBlocks_];
            for(register label i=0;i<numBlocks_;++i)
                dataptr1[i] = dataPtr_[i];
            
            if( dataPtr_ )
                delete [] dataPtr_;
            dataPtr_ = dataptr1;
        }
        
        for(register label i=numBlocks_;i<numblock1;++i)
            dataPtr_[i] = new T[blockSize];
    }
    
    numBlocks_ = numblock1;
    N_ = numBlocks_ * blockSize;
}

template<class T, Foam::label Offset>
void Foam::LongList<T, Offset>::clearOut()
{
    for(register label i=0;i<numBlocks_;++i)
        delete [] dataPtr_[i];
    
    if( dataPtr_ )
    {
        delete [] dataPtr_;
        dataPtr_ = NULL;
    }
    
    N_ = 0;
    numBlocks_ = 0;
    numAllocatedBlocks_ = 0;
    nextFree_ = 0;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template<class T, Foam::label Offset>
inline Foam::LongList<T, Offset>::LongList()
:
    N_(0),
    nextFree_(0),
    numBlocks_(0),
    numAllocatedBlocks_(0),
    shift_(),
    mask_(),
    dataPtr_(NULL)
{
    initializeParameters();
}

//- Construct given size
template<class T, Foam::label Offset>
inline Foam::LongList<T, Offset>::LongList(const label s)
:
    N_(0),
    nextFree_(0),
    numBlocks_(0),
    numAllocatedBlocks_(0),
    shift_(),
    mask_(),
    dataPtr_(NULL)
{
    initializeParameters();
    setSize(s);
}


//- Construct given size
template<class T, Foam::label Offset>
inline Foam::LongList<T, Offset>::LongList(const label s, const T& t)
:
    N_(0),
    nextFree_(0),
    numBlocks_(0),
    numAllocatedBlocks_(0),
    shift_(),
    mask_(),
    dataPtr_(NULL)
{
    initializeParameters();
    setSize(s);
    *this = t;
}

template<class T, Foam::label Offset>
inline Foam::LongList<T, Offset>::LongList(const LongList<T, Offset>& ol)
:
    N_(0),
    nextFree_(0),
    numBlocks_(0),
    numAllocatedBlocks_(0),
    shift_(ol.shift_),
    mask_(ol.mask_),
    dataPtr_(NULL)
{
    *this = ol;
}

template<class T, Foam::label Offset>
inline Foam::LongList<T, Offset>::~LongList()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, Foam::label Offset>
inline Foam::label Foam::LongList<T, Offset>::size() const
{
    return nextFree_;
}

template<class T, Foam::label Offset>
inline Foam::label Foam::LongList<T, Offset>::byteSize() const
{
    if( !contiguous<T>() )
    {
        FatalErrorIn("LongList<T, Offset>::byteSize()")
            << "Cannot return the binary size of a list of "
               "non-primitive elements"
            << abort(FatalError);
    }

    return nextFree_*sizeof(T);
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::setSize(const label i)
{
    allocateSize(i);
    nextFree_ = i;
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::clear()
{
    nextFree_ = 0;
}


template<class T, Foam::label Offset>
inline Foam::LongList<T, Offset>&
Foam::LongList<T, Offset>::shrink()
{
    setSize(nextFree_);
    return *this;
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::transfer(LongList<T, Offset>& ol)
{
    clearOut();
    dataPtr_ = ol.dataPtr_;
    N_ = ol.N_;
    nextFree_ = ol.nextFree_;
    numBlocks_ = ol.numBlocks_;
    numAllocatedBlocks_ = ol.numAllocatedBlocks_;
    shift_ = ol.shift_;
    mask_ = ol.mask_;
    
    ol.dataPtr_ = NULL;
    ol.N_ = 0;
    ol.nextFree_ = 0;
    ol.numBlocks_ = 0;
    ol.numAllocatedBlocks_ = 0;
}


template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::append(const T& e)
{
    if( nextFree_ >= N_ )
    {
        allocateSize(nextFree_+1);
    }
    
    operator[](nextFree_++) = e;
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::appendIfNotIn(const T& e)
{
    if( !contains(e) )
         append(e);
}

template<class T, Foam::label Offset>
inline bool Foam::LongList<T, Offset>::contains(const T& e) const
{
    for(register label i=0;i<nextFree_;++i)
        if( (*this)[i] == e )
            return true;
        
    return false;
}

template<class T, Foam::label Offset>
inline Foam::label Foam::LongList<T, Offset>::containsAtPosition
(
    const T& e
) const
{
    for(register label i=0;i<nextFree_;++i)
        if( (*this)[i] == e )
            return i;
        
    return -1;
}

template<class T, Foam::label Offset>
inline T Foam::LongList<T, Offset>::remove(const label i)
{
    if( nextFree_ == 0 )
    {
        FatalErrorIn
        (
            "void Foam::LongList<T, Offset>::remove()"
        )   << "List is empty" << abort(FatalError);
    }
    
    T el = operator[](i);
    operator[](i) = operator[](nextFree_-1);
    --nextFree_;
    return el;
}

template<class T, Foam::label Offset>
inline T Foam::LongList<T, Offset>::removeLastElement()
{
    if( nextFree_ == 0 )
    {
        FatalErrorIn
        (
            "void Foam::LongList<T, Offset>::remove()"
        )   << "List is empty" << abort(FatalError);
    }
    
    T lastEl = operator[](nextFree_-1);
    --nextFree_;
    return lastEl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, Foam::label Offset>
inline const T& Foam::LongList<T, Offset>::operator[](const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
        
    return dataPtr_[i>>shift_][i&mask_];
}

template<class T, Foam::label Offset>
inline T& Foam::LongList<T, Offset>::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
        
    return dataPtr_[i>>shift_][i&mask_];
}

template<class T, Foam::label Offset>
inline T& Foam::LongList<T, Offset>::operator()(const label i)
{
    if( i >= nextFree_ )
        setSize(i+1);

    return operator[](i);
}


template<class T, Foam::label Offset>
inline T& Foam::LongList<T, Offset>::newElmt(const label i)
{
    return operator()(i);
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::operator=(const T& t)
{
    for(register label i=0;i<nextFree_;++i)
        operator[](i) = t;
}

template<class T, Foam::label Offset>
inline void Foam::LongList<T, Offset>::operator=(const LongList<T, Offset>& l)
{
    setSize(l.size());
    for(register label i=0;i<l.nextFree_;++i)
        operator[](i) = l[i];
}


// ************************************************************************* //
