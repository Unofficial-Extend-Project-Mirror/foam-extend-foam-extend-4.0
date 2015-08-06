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

template<class graphType>
inline void Foam::graphRow<graphType>::checkIndex(const label i) const
{
    if( (i < 0) || (i >=data_.sizeOfRow(rowI_)) )
    {
        FatalErrorIn
        (
            "void Foam::graphRow<graphType>::"
            "checkIndex(const label i) const"
        ) << "Row index " << rowI_
            << " is not in range " << Foam::label(0)
            << " and " << data_.sizeOfRow(rowI_) << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class graphType>
inline Foam::graphRow<graphType>::graphRow(graphType& g, const label i)
:
    data_(g),
    rowI_(i)
{
}

template<class graphType>
inline Foam::graphRow<graphType>::graphRow
(
    const graphRow<graphType>& ol
)
:
    data_(ol.data_),
    rowI_(ol.rowI_)
{
}

template<class graphType>
inline Foam::graphRow<graphType>::~graphRow()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class graphType>
inline Foam::label Foam::graphRow<graphType>::size() const
{
    return data_.sizeOfRow(rowI_);
}

template<class graphType>
inline void Foam::graphRow<graphType>::setSize(const label s)
{
    data_.setRowSize(rowI_, s);
}

template<class graphType>
inline void Foam::graphRow<graphType>::clear()
{
    data_.setRowSize(rowI_, 0);
}

template<class graphType>
inline void Foam::graphRow<graphType>::append(const label el)
{
    data_.append(rowI_, el);
}

template<class graphType>
inline void Foam::graphRow<graphType>::appendIfNotIn(const label el)
{
    data_.appendIfNotIn(rowI_, el);
}

template<class graphType>
inline bool Foam::graphRow<graphType>::contains(const label e) const
{
    return data_.contains(rowI_, e);
}

template<class graphType>
inline Foam::label Foam::graphRow<graphType>::containsAtPosition
(
    const label e
) const
{
    return data_.containsAtPosition(rowI_, e);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class graphType>
inline Foam::label Foam::graphRow<graphType>::operator[](const label i) const
{
    return data_(rowI_, i);
}

template<class graphType>
inline Foam::label& Foam::graphRow<graphType>::operator[](const label i)
{
    return data_(rowI_, i);
}

template<class graphType>
inline void Foam::graphRow<graphType>::operator=
(
    const graphRow<graphType>& l
)
{
    data_.setRowSize(rowI_, l.size());
    for(label i=0;i<l.size();++i)
        data_(rowI_, i) = l[i];
}

template<class graphType>
template<class listType>
inline void Foam::graphRow<graphType>::operator=(const listType& l)
{
    data_.setRowSize(rowI_, l.size());
    for(label i=0;i<l.size();++i)
        data_(rowI_, i) = l[i];
}

template<class graphType>
inline Foam::Ostream& operator<<
(
    Foam::Ostream& os,
    const Foam::graphRow<graphType>& r
)
{
    os << r.size() << "(";
    for(Foam::label i=0;i<r.size();++i)
        os << r[i] << " ";
    os << ")";
    
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
