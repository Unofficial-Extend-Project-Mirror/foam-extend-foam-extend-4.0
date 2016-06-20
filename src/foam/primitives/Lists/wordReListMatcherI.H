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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wordReListMatcher::wordReListMatcher
(
    const UList<wordRe>& lst
)
:
    reList_(lst)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::wordReListMatcher::size() const
{
    return reList_.size();
}


inline bool Foam::wordReListMatcher::empty() const
{
    return reList_.empty();
}


inline const Foam::UList<Foam::wordRe>&
Foam::wordReListMatcher::operator()() const
{
    return reList_;
}


inline bool Foam::wordReListMatcher::match
(
    const string& str,
    bool literalMatch
) const
{
    const label nElem = reList_.size();
    for (label elemI = 0; elemI < nElem; ++elemI)
    {
        if (reList_[elemI].match(str, literalMatch))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
