/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "profilingInfo.H"

#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::profilingInfo::nextId_(0);

Foam::label Foam::profilingInfo::getID()
{
    nextId_++;
    return nextId_;
}

void Foam::profilingInfo::raiseID(label maxVal)
{
    if(maxVal>nextId_) {
        nextId_=maxVal;
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingInfo::profilingInfo()
:
    calls_(0),
    totalTime_(0.),
    childTime_(0.),
    id_(getID()),
    parent_(*this),
    description_("application::main"),
    onStack_(false)
{}


Foam::profilingInfo::profilingInfo(profilingInfo &parent,const string &descr)
:
    calls_(0),
    totalTime_(0.),
    childTime_(0.),
    id_(getID()),
    parent_(parent),
    description_(descr),
    onStack_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingInfo::~profilingInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingInfo::update(scalar elapsedTimee)
{
    calls_++;
    totalTime_+=elapsedTimee;
    if(id()!=parent().id()) {
        parent_.childTime_+=elapsedTimee;
    }
}

void Foam::profilingInfo::writeWithOffset(Ostream &os,bool offset,scalar time,scalar childTimes) const
{
    dictionary tmp;
    
    tmp.add("id",id());
    if(id()!=parent().id()) {
        tmp.add("parentId",parent().id());
    }
    tmp.add("description",description());
    tmp.add("calls",calls()+(offset ? 1 : 0));
    tmp.add("totalTime",totalTime()+time);
    tmp.add("childTime",childTime()+childTimes);
    tmp.add("onStack",onStack());

    os << tmp;  
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const profilingInfo& info)
{
    info.writeWithOffset(os);

    return os;
}

// ************************************************************************* //
