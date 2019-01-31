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

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void Foam::VRWGraphList::checkIndex
(
    const label i,
    const label j,
    const label k
) const
{
    if( (i < 0) || (i >= rows_.size()) )
    {
        FatalErrorIn
        (
            "void Foam::VRWGraphList<T,width>::"
            "checkIndex(const label i, const label j, const label k) const"
        ) << "Position index " << Foam::label(i)
            << " is not in range " << Foam::label(0)
            << " and " << rows_.size() << abort(FatalError);
    }
    
    if( (j < 0) || (j >= rows_[i].size()) )
        FatalErrorIn
        (
            "void Foam::VRWGraphList<T,width>::"
            "checkIndex(label const i, const label j, const label k) const"
        ) << "Row index " << Foam::label(j)
            << " is not in range " << Foam::label(0)
            << " and " << rows_[i].size() << abort(FatalError);
    
    if( (k < 0) || (k >= data_.sizeOfRow(rows_[i].start()+j)) )
        FatalErrorIn
        (
            "void Foam::VRWGraphList<T,width>::"
            "checkIndex(label const i, const label j, const label k) const"
        ) << "Data index " << Foam::label(k)
            << " is not in range " << Foam::label(0)
            << " and " << data_.sizeOfRow(rows_[i].start()+j)
            << abort(FatalError);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline Foam::VRWGraphList::VRWGraphList()
:
    data_(),
    rows_()
{
}

inline Foam::VRWGraphList::VRWGraphList
(
    const VRWGraphList& ol
)
:
    data_(ol.data_),
    rows_(ol.rows_)
{
}

inline Foam::VRWGraphList::~VRWGraphList()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::VRWGraphList::size() const
{
    return rows_.size();
}

inline Foam::label Foam::VRWGraphList::sizeOfGraph(const label posI) const
{
    return rows_[posI].size();
}

inline Foam::label Foam::VRWGraphList::sizeOfRow
(
    const label posI,
    const label rowI
) const
{
    return data_.sizeOfRow(rows_[posI].start()+rowI);
}

inline void Foam::VRWGraphList::clear()
{
    data_.setSize(0);
    rows_.setSize(0);
}

template<class GraphType>
inline void Foam::VRWGraphList::appendGraph
(
    const GraphType& l
)
{
    rowElement re(data_.size(), l.size());
    
    for(label i=0;i<l.size();++i)
        data_.appendList(l[i]);
    
    rows_.append(re);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::VRWGraphList::operator()
(
    const label i,
    const label j,
    const label k
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j, k);
    #endif
        
    return data_(rows_[i].start() + j, k);
}


inline Foam::label& Foam::VRWGraphList::operator()
(
    const label i, const label j, const label k
)
{
    #ifdef FULLDEBUG
    checkIndex(i, j, k);
    #endif
        
    return data_(rows_[i].start() + j, k);
}

inline const Foam::subGraph<const VRWGraph> Foam::VRWGraphList::operator[]
(
    const label i
) const
{
    return Foam::subGraph<const VRWGraph>(data_, i, rows_[i].size());
}

inline void Foam::VRWGraphList::operator=(const VRWGraphList& l)
{
    data_ = l.data_;
    rows_ = l.rows_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
