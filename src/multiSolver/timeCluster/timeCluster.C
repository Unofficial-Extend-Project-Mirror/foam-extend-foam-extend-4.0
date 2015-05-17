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

#include "timeClusterList.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::timeCluster::typeName = "timeCluster";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeCluster::timeCluster()
{}


Foam::timeCluster::timeCluster
(
    const instantList& times,
    const scalar globalOffset,
    const label globalIndex,
    const label superLoop,
    const word& solverDomainName,
    const word& preConName
)
:
    instantList(times),
    globalOffset_(globalOffset),
    globalIndex_(globalIndex),
    superLoop_(superLoop),
    solverDomainName_(solverDomainName),
    preConName_(preConName)
{}

Foam::timeCluster::timeCluster
(
    const timeCluster& tc,
    const label index
)
:
    instantList(1, tc[index]),
    globalOffset_(tc.globalOffset_),
    globalIndex_(tc.globalIndex_),
    superLoop_(tc.superLoop_),
    solverDomainName_(tc.solverDomainName_),
    preConName_(tc.preConName_)
{}


Foam::timeCluster::timeCluster(const Foam::scalar t)
:
    instantList(1, instant(0)),
    globalOffset_(0),
    globalIndex_(0),
    superLoop_(0),
    solverDomainName_(word::null),
    preConName_(word::null)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::timeCluster::globalValue(const label& index) const
{
    return this->operator[](index).value() + globalOffset_;
}

Foam::scalar Foam::timeCluster::globalMinValue() const
{
    return this->localMinValue() + globalOffset_;
}


Foam::scalar Foam::timeCluster::globalMaxValue() const
{
    return this->localMaxValue() + globalOffset_;
}


Foam::label Foam::timeCluster::globalMinIndex() const
{
    return this->localMinIndex();
}


Foam::label Foam::timeCluster::globalMaxIndex() const
{
    return this->localMaxIndex();
}


Foam::scalar Foam::timeCluster::globalFindClosestTimeValue
(
    const scalar timeValue
) const
{
    return this->operator[]
    (
        Foam::Time::findClosestTimeIndex
        (
            this->times(), timeValue - globalOffset_
        )
    ).value() + globalOffset_;
}


Foam::label Foam::timeCluster::globalFindClosestTimeIndex
(
    const scalar timeValue
) const
{
    return Foam::Time::findClosestTimeIndex
    (
        this->times(), timeValue - globalOffset_
    );
}


Foam::scalar Foam::timeCluster::localValue(const label& index) const
{
    return this->operator[](index).value();
}


Foam::scalar Foam::timeCluster::localMinValue() const
{
    return this->operator[](this->localMinIndex()).value();
}


Foam::scalar Foam::timeCluster::localMaxValue() const
{
    return this->operator[](this->localMaxIndex()).value();
}


Foam::label Foam::timeCluster::localMinIndex() const
{
    label bestIndex(0);
    scalar min(VGREAT);
    forAll(*this, i)
    {
        if (this->operator[](i).value() < min)
        {
            min = this->operator[](i).value();
            bestIndex = i;
        }
    }
    return bestIndex;
}


Foam::label Foam::timeCluster::localMaxIndex() const
{
    label bestIndex(0);
    scalar max(0);
    forAll(*this, i)
    {
        if (this->operator[](i).value() > max)
        {
            max = this->operator[](i).value();
            bestIndex = i;
        }
    }
    return bestIndex;
}


Foam::scalar Foam::timeCluster::localFindClosestTimeValue
(
    const scalar timeValue
) const
{
    return this->operator[]
    (
        Foam::Time::findClosestTimeIndex
        (
            this->times(), timeValue
        )
    ).value();
}


Foam::label Foam::timeCluster::localFindClosestTimeIndex
(
    const scalar timeValue
) const
{
    return Foam::Time::findClosestTimeIndex
    (
        this->times(), timeValue
    );
}


Foam::timeCluster Foam::timeCluster::operator()(const Foam::label index) const
{
    return timeCluster(*this, index);
}


// * * * * * * * * * * * * * Friend IOstream Operators * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, timeCluster& I)
{
    return is   >> I.globalOffset_
                >> I.globalIndex_
                >> I.superLoop_
                >> I.solverDomainName_
                >> I.preConName_
                >> I.times();
}


Foam::Ostream& Foam::operator<<(Ostream& os, const timeCluster& I)
{
    return os   << "/* globalOffset: */\t" << I.globalOffset_ << nl
                << "/* globalIndex:  */\t" << I.globalIndex_ << nl
                << "/* superLoop:    */\t" << I.superLoop_ << nl
                << "/* solverDomain: */\t" << I.solverDomainName_ << nl
                << "/* preConName:   */\t" << I.preConName_ << nl
                << "/* Instant list: */\t" << I.times();
}

// ************************************************************************* //
