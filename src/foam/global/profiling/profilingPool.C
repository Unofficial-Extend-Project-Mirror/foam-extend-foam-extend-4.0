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

#include "profilingPool.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::profilingPool* Foam::profilingPool::thePool_(NULL);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingPool::profilingPool(
    const IOobject &ob,
    const Time &owner
)
    :
    regIOobject(ob),
    globalTime_(),
    owner_(owner)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingPool::~profilingPool()
{
    for(mapIterator it = map().begin(); it != map().end(); ++it)
    {
        delete it->second;
    }

    map().erase(allInfo_.begin(), allInfo_.end());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingPool::initProfiling(
    const IOobject &ob,
    const Time &owner
)
{
    if (!thePool_)
    {
        thePool_ = new profilingPool(ob,owner);
        profilingInfo *master=new profilingInfo();
        thePool_->map().insert(make_pair(master->description(),master));
        thePool_->stack().push(*master);
        profilingPool::rememberTimer(*master,thePool_->globalTime_);
    }
}

void Foam::profilingPool::stopProfiling(
    const Time &owner
)
{
    if (thePool_ && (&owner)==&(thePool_->owner()))
    {
        delete thePool_;
        thePool_=NULL;
    }
}

Foam::profilingInfo &Foam::profilingPool::getInfo(const string& name)
{
    if (!thePool_)
    {
        FatalErrorIn("profilingPool::addInfo(const string& name)")
            << "Singleton not initialized\n" << endl
            << abort(FatalError);
    }

    profilingStack& stack = thePool_->stack();
    mapType& map = thePool_->map();

    profilingInfo* found = NULL;

    for
    (
        mapIterator it = map.lower_bound(name);
        it != map.upper_bound(name);
        ++it
    )
    {
        if (it->second->parent().id()==stack.top().id())
        {
            found = it->second;
            break;
        }
    }

    if (!found)
    {
        found = new profilingInfo(stack.top(),name);

        map.insert(make_pair(name,found));
    }

    stack.push(*found);
    return *found;
}


void Foam::profilingPool::rememberTimer
(
    const profilingInfo& info,
    clockTime& timer
)
{
    if(!thePool_)
    {
        FatalErrorIn
        (
            "profilingPool::rememberTimer(const profilingInfo Foam&info, "
            "clockTime& timer)"
        )   << "Singleton not initialized\n" << endl
            << abort(FatalError);
    }

    thePool_->stack().addTimer(info, timer);
}


void Foam::profilingPool::remove(const profilingInfo &info)
{
    if(!thePool_)
    {
        FatalErrorIn("profilingPool::addInfo(const string& name)")
            << "Singleton not initialized\n" << endl
            << abort(FatalError);
    }

    profilingStack& stack = thePool_->stack();

    if(info.id() != stack.top().id())
    {
        FatalErrorIn("profilingPool::update(const string &name)")
            << "The id " << info.id() << " of the updated info "
            << info.description()
            << " is no the same as the one on top of the stack: "
            << stack.top().id() << " (" << stack.top().description()
            << ")\n" << endl
            << abort(FatalError);
    }

    stack.pop();
}


bool Foam::profilingPool::writeData(Ostream& os) const
{
    os  << "profilingInfo" << nl << indent
        << token::BEGIN_LIST << incrIndent << nl;

    stack().writeStackContents(os);

    for(mapConstIterator it = map().begin(); it != map().end(); ++it)
    {
        if(!it->second->onStack())
        {
            os << *(it->second);
        }
    }

    os  << decrIndent << indent << token::END_LIST
        << token::END_STATEMENT << endl;

    return os;
}

// ************************************************************************* //
