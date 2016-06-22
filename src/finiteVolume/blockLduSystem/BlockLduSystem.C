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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const lduMesh& ldu
)
:
    BlockLduMatrix<blockType>(ldu),
    source_(ldu.lduAddr().size(), pTraits<sourceType>::zero)
{}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const lduMesh& ldu,
    const Field<sourceType>& s
)
:
    BlockLduMatrix<blockType>(ldu),
    source_(s)
{
    if (ldu.lduAddr().size() != s.size())
    {
        FatalErrorIn
        (
            "BlockLduSystem::BlockLduSystem\n"
            "(\n"
            "    const lduMesh& ldu,"
            "    const Field<sourceType>& s,"
            ")\n"
        )   << "Sizes of ldu addressing and source field are not the same."
            << abort(FatalError);
    }
}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const BlockLduMatrix<blockType>& bm,
    const Field<sourceType>& s
)
:
    BlockLduMatrix<blockType>(bm),
    source_(s)
{
    if (this->lduAddr().size() != s.size())
    {
        FatalErrorIn
        (
            "BlockLduSystem::BlockLduSystem\n"
            "(\n"
            "    const BlockLduMatrix<blockType>& bm,"
            "    const Field<sourceType>& s,"
            ")\n"
        )   << "Sizes of block matrix and source field are not the same."
            << abort(FatalError);
    }
}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const BlockLduSystem<blockType, sourceType>& bs
)
:
    BlockLduMatrix<blockType>(bs),
    source_(bs.source())
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class blockType, class sourceType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    os  << static_cast<const BlockLduMatrix<blockType>&>(bs) << nl
        << bs.source() << endl;

    return os;
}


// ************************************************************************* //
