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

#include "ProfilingPool.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::ProfilingPool* Foam::ProfilingPool::thePool_(NULL);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ProfilingPool::ProfilingPool(const IOobject &ob)
    :
    regIOobject(ob),
    globalTime_()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ProfilingPool::~ProfilingPool()
{
    for(mapIterator it=map().begin();it!=map().end();++it) {
        delete it->second;
    }
    map().erase(allInfo_.begin(),allInfo_.end());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ProfilingPool::initProfiling(const IOobject &ob)
{
    if(thePool_!=NULL) {
        WarningIn("Foam::ProfilingPool::initProfiling(const IOobject &)")
            << "Singleton already initialized\n" << endl;
    } else {
        thePool_=new ProfilingPool(ob);
        ProfilingInfo *master=new ProfilingInfo();
        thePool_->map().insert(make_pair(master->description(),master));
        thePool_->stack().push(*master);
        ProfilingPool::rememberTimer(*master,thePool_->globalTime_);
    }
}

Foam::ProfilingInfo &Foam::ProfilingPool::getInfo(const string &name)
{
    if(thePool_==NULL) {
        FatalErrorIn("Foam::ProfilingPool::addInfo(const string &name)")
            << "Sinleton not initialized\n" << endl
                << abort(FatalError);
    } 

    ProfilingStack &stack=thePool_->stack();
    mapType &map=thePool_->map();

    ProfilingInfo *found=NULL;

    for(mapIterator it=map.lower_bound(name);it!=map.upper_bound(name);++it) {
        if(it->second->parent().id()==stack.top().id()) {
            found=it->second;
            break;
        }
    }

    if(found==NULL) {
        found=new ProfilingInfo(stack.top(),name);

        map.insert(make_pair(name,found));
    }

    stack.push(*found);
    return *found;
}

void Foam::ProfilingPool::rememberTimer(const ProfilingInfo &info,clockTime &timer)
{
    if(thePool_==NULL) {
        FatalErrorIn("Foam::ProfilingPool::rememberTimer(const ProfilingInfo &info,clockTime &timer)")
            << "Singleton not initialized\n" << endl
                << abort(FatalError);
    } 
    
    thePool_->stack().addTimer(info,timer);
}

void Foam::ProfilingPool::remove(const ProfilingInfo &info)
{
    if(thePool_==NULL) {
        FatalErrorIn("Foam::ProfilingPool::addInfo(const string &name)")
            << "Singleton not initialized\n" << endl
                << abort(FatalError);
    } 

    ProfilingStack &stack=thePool_->stack();

    if(info.id()!=stack.top().id()) {
        FatalErrorIn("Foam::ProfilingPool::update(const string &name)")
            << "The id " << info.id() << " of the updated info " << info.description()
                << " is no the same as the one on top of the stack: " 
                << stack.top().id() << " (" << stack.top().description() << ")\n" << endl
                << abort(FatalError);
    }

    stack.pop();
}

bool Foam::ProfilingPool::writeData(Ostream &os) const
{
    os << "profilingInfo" << nl << indent << token::BEGIN_LIST << incrIndent << nl;

    stack().writeStackContents(os);

    for(mapConstIterator it=map().begin();it!=map().end();++it) {
        if(!it->second->onStack()) {
            os << *(it->second);
        }
    }

    os << decrIndent << indent << token::END_LIST << token::END_STATEMENT << endl;

    return os;
}

// ************************************************************************* //
