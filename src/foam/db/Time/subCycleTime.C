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

#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subCycleTime::subCycleTime(Time& t, const label nSubCycles)
:
    time_(t),
    nSubCycles_(nSubCycles),
    subCycleIndex_(0)
{
    time_.subCycle(nSubCycles_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subCycleTime::~subCycleTime()
{
    endSubCycle();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::subCycleTime::end() const
{
    return subCycleIndex_ > nSubCycles_;
}


void Foam::subCycleTime::endSubCycle()
{
    time_.endSubCycle();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::subCycleTime& Foam::subCycleTime::operator++()
{
    time_++;
    subCycleIndex_++;
    return *this;
}


Foam::subCycleTime& Foam::subCycleTime::operator++(int)
{
    return operator++();
}


// ************************************************************************* //
