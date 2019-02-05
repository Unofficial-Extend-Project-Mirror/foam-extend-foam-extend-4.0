/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

namespace Foam
{

template<class graphType>
inline void Foam::subGraph<graphType>::checkIndex(const label i) const
{
    if( (i < 0) || (i >= size_) )
    {
        FatalErrorIn
        (
            "void Foam::subGraph<graphType>::"
            "checkIndex(const label i) const"
        ) << "Row index " << i
            << " is not in range " << Foam::label(0)
            << " and " << size_ << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class graphType>
inline Foam::subGraph<graphType>::subGraph
(
    graphType& g,
    const label start,
    const label size
)
:
    data_(g),
    start_(start),
    size_(size)
{
}

template<class graphType>
inline Foam::subGraph<graphType>::subGraph
(
    const subGraph<graphType>& ol
)
:
    data_(ol.data_),
    start_(ol.start_),
    size_(ol.size_)
{
}

template<class graphType>
inline Foam::subGraph<graphType>::~subGraph()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class graphType>
inline Foam::label Foam::subGraph<graphType>::size() const
{
    return size_;
}

template<class graphType>
inline Foam::label Foam::subGraph<graphType>::sizeOfRow(const label rowI) const
{
    return data_.sizeOfRow(start_+rowI);
}

template<class graphType>
inline void Foam::subGraph<graphType>::append(const label rowI, const label el)
{
    data_.append(start_+rowI, el);
}

template<class graphType>
inline void Foam::subGraph<graphType>::appendIfNotIn
(
    const label rowI,
    const label el
)
{
    data_.appendIfNotIn(start_+rowI, el);
}

template<class graphType>
inline bool Foam::subGraph<graphType>::contains
(
    const label rowI,
    const label e
) const
{
    return data_.contains(start_+rowI, e);
}

template<class graphType>
inline Foam::label Foam::subGraph<graphType>::containsAtPosition
(
    const label rowI,
    const label e
) const
{
    return data_.containsAtPosition(start_+rowI, e);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class graphType>
inline Foam::label Foam::subGraph<graphType>::operator()
(
    const label i,
    const label j
) const
{
    # ifdef FULLDEBUG
    checkIndex(i);
    # endif

    return data_(start_+i, j);
}

template<class graphType>
inline Foam::label& Foam::subGraph<graphType>::operator()
(
    const label i,
    const label j
)
{
    return data_(start_+i, j);
}


template<class graphType>
inline const Foam::graphRow<const graphType>
Foam::subGraph<graphType>::operator[]
(
    const label i
) const
{
    return data_[start_+i];
}

template<class graphType>
inline Foam::graphRow<graphType>
Foam::subGraph<graphType>::operator[](const label i)
{
    return data_[start_+i];
}

template<class graphType>
inline Foam::Ostream& operator<<
(
    Foam::Ostream& os,
    const Foam::subGraph<graphType>& sg
)
{
    os << sg.size() << "\n" << "(";
    for(Foam::label i=0;i<sg.size();++i)
    {
        os << "\n" << sg.sizeOfRow(i) << "(";
        for(Foam::label j=0;j<sg.sizeOfRow(i);++j)
        {
            if( j > 0 ) os << " ";

            os << sg(i, j);
        }

        os << ")";
    }

    os << "\n" << ")";

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
