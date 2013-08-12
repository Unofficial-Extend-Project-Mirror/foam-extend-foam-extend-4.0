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

Class
    Foam::ProfilingPool

Description
    Collects all the data for Profiling 

SourceFiles
    ProfilingPool.C

\*---------------------------------------------------------------------------*/

#ifndef ProfilingPool_H
#define ProfilingPool_H

#include "regIOobject.H"
#include "clockTime.H"

#include <map>

#include "ProfilingInfo.H"
#include "ProfilingStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
// class Istream;
class Ostream;

// // Forward declaration of friend functions and operators
// Istream& operator>>(Istream&, ProfilingPool&);
// Ostream& operator<<(Ostream&, const ProfilingPool&);


/*---------------------------------------------------------------------------*\
                           Class ProfilingPool Declaration
\*---------------------------------------------------------------------------*/

class ProfilingPool
:
    public regIOobject
{
    // Private data

    typedef std::multimap<Foam::string,Foam::ProfilingInfo*> mapType;
    typedef std::pair<Foam::string,Foam::ProfilingInfo*> mapValues;
    typedef mapType::iterator mapIterator;
    typedef mapType::const_iterator mapConstIterator;

    mapType allInfo_;

    ProfilingStack theStack_;

    clockTime globalTime_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    ProfilingPool(const ProfilingPool&);
    
    //- Disallow default bitwise assignment
    void operator=(const ProfilingPool&);
    
    // Static data members

    //- the only possible Pool-Object
    static ProfilingPool *thePool_;

    //- Construct null
    ProfilingPool(const IOobject &);

    // Destructor
    
    ~ProfilingPool();

    mapType &map() 
        { return allInfo_; }

    const mapType &map() const
        { return allInfo_; }

    ProfilingStack &stack() 
        { return theStack_; }

    const ProfilingStack &stack() const
        { return theStack_; }

public:

    static void initProfiling(const IOobject &);

    static ProfilingInfo &getInfo(const string &name);

    static void remove(const ProfilingInfo &info);

    static void rememberTimer(const ProfilingInfo &info,clockTime &timer);

    virtual bool writeData(Ostream &) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
