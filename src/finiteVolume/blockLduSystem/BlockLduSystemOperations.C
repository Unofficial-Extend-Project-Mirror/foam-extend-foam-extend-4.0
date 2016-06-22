/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlockLduSystem.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::negate()
{
    BlockLduMatrix<blockType>::negate();
    source_.negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator=
(
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    if (this == &bs)
    {
        FatalErrorIn
        (
            "void BlockLduSystem<blockType, sourceType>::operator="
            "(const BlockLduSystem<blockType, sourceType>& bs)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    BlockLduMatrix<blockType>::operator=(bs);
    source_ = bs.source();
}


template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator+=
(
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    BlockLduMatrix<blockType>::operator+=(bs);
    source_ += bs.source();
}

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator-=
(
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    BlockLduMatrix<blockType>::operator-=(bs);
    source_ -= bs.source();
}

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator*=
(
    const scalarField& sf
)
{
    BlockLduMatrix<blockType>::operator*=(sf);
    source_ *= sf;
}

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator*=
(
    const scalar s
)
{
    BlockLduMatrix<blockType>::operator*=(s);
    source_ *= s;
}


// ************************************************************************* //
