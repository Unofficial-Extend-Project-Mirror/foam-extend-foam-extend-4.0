/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const instantList& timeCluster::times() const
{
    const instantList& times(*this);
    return times;
}

inline instantList& timeCluster::times()
{
    instantList& times(*this);
    return times;
}

inline scalar timeCluster::globalOffset() const
{
    return globalOffset_;
}

inline scalar& timeCluster::globalOffset()
{
    return globalOffset_;
}

inline label timeCluster::globalIndex() const
{
    return globalIndex_;
}

inline label& timeCluster::globalIndex()
{
    return globalIndex_;
}

inline label timeCluster::superLoop() const
{
    return superLoop_;
}

inline label& timeCluster::superLoop()
{
    return superLoop_;
}

inline const word& timeCluster::solverDomainName() const
{
    return solverDomainName_;
}

inline word& timeCluster::solverDomainName()
{
    return solverDomainName_;
}

inline const word& timeCluster::preConName() const
{
    return preConName_;
}

inline word& timeCluster::preConName()
{
    return preConName_;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
