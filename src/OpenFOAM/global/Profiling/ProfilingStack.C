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

#include "ProfilingStack.H"
#include "ProfilingInfo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ProfilingStack::ProfilingStack()
:
    LIFOStack<ProfilingInfo*>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ProfilingStack::~ProfilingStack()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::ProfilingInfo &Foam::ProfilingStack::top() const
{
    return *LIFOStack<ProfilingInfo*>::top();
}

Foam::ProfilingInfo &Foam::ProfilingStack::bottom() const
{
    return *LIFOStack<ProfilingInfo*>::bottom();
}

bool Foam::ProfilingStack::empty() const
{
    return LIFOStack<ProfilingInfo*>::empty();
}

void Foam::ProfilingStack::push(ProfilingInfo &a)
{
    LIFOStack<ProfilingInfo*>::push(&a);
    top().addedToStack();
}

Foam::ProfilingInfo &Foam::ProfilingStack::pop()
{
    top().removedFromStack();
    return *LIFOStack<ProfilingInfo*>::pop();
}

void Foam::ProfilingStack::writeStackContents(Ostream &os) const
{
    if(empty()) {
        return;
    }
    const_iterator it=begin();
    scalar oldElapsed=0;
    do {
        const ProfilingInfo &info=*(*it);
        scalar elapsed=timers_[info.id()]->elapsedTime();

        info.writeWithOffset(os,true,elapsed,oldElapsed);

        oldElapsed=elapsed;
        ++it;
    } while(it!=end());
}

void Foam::ProfilingStack::addTimer(const ProfilingInfo &info,clockTime &timer)
{
    timers_.insert(info.id(),&timer);
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
