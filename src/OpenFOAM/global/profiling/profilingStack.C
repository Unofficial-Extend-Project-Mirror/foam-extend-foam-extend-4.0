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

#include "profilingStack.H"
#include "profilingInfo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingStack::profilingStack()
:
    LIFOStack<profilingInfo*>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingStack::~profilingStack()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::profilingInfo &Foam::profilingStack::top() const
{
    return *LIFOStack<profilingInfo*>::top();
}

Foam::profilingInfo &Foam::profilingStack::bottom() const
{
    return *LIFOStack<profilingInfo*>::bottom();
}

bool Foam::profilingStack::empty() const
{
    return LIFOStack<profilingInfo*>::empty();
}

void Foam::profilingStack::push(profilingInfo &a)
{
    LIFOStack<profilingInfo*>::push(&a);
    top().addedToStack();
}

Foam::profilingInfo &Foam::profilingStack::pop()
{
    top().removedFromStack();
    return *LIFOStack<profilingInfo*>::pop();
}

void Foam::profilingStack::writeStackContents(Ostream &os) const
{
    if(empty()) {
        return;
    }
    const_iterator it=begin();
    scalar oldElapsed=0;
    do {
        const profilingInfo &info=*(*it);
        scalar elapsed=timers_[info.id()]->elapsedTime();

        info.writeWithOffset(os,true,elapsed,oldElapsed);

        oldElapsed=elapsed;
        ++it;
    } while(it!=end());
}

void Foam::profilingStack::addTimer(const profilingInfo &info,clockTime &timer)
{
    timers_.insert(info.id(),&timer);
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
