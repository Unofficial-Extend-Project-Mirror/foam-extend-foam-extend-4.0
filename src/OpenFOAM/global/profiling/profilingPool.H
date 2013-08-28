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
    Foam::profilingPool

Description
    Collects all the data for profiling 

SourceFiles
    profilingPool.C

\*---------------------------------------------------------------------------*/

#ifndef profilingPool_H
#define profilingPool_H

#include "regIOobject.H"
#include "clockTime.H"

#include <map>

#include "profilingInfo.H"
#include "profilingStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
// class Istream;
class Ostream;

// // Forward declaration of friend functions and operators
// Istream& operator>>(Istream&, profilingPool&);
// Ostream& operator<<(Ostream&, const profilingPool&);


/*---------------------------------------------------------------------------*\
                           Class profilingPool Declaration
\*---------------------------------------------------------------------------*/

class profilingPool
:
    public regIOobject
{
    // Private typedefs

        typedef std::multimap<Foam::string,Foam::profilingInfo*> mapType;
        typedef std::pair<Foam::string,Foam::profilingInfo*> mapValues;
        typedef mapType::iterator mapIterator;
        typedef mapType::const_iterator mapConstIterator;


    // Private data

        mapType allInfo_;

        profilingStack theStack_;

        clockTime globalTime_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        profilingPool(const profilingPool&);
    
        //- Disallow default bitwise assignment
        void operator=(const profilingPool&);
    

    // Static data members

        //- The only possible Pool-Object
        static profilingPool *thePool_;


    // Constructors

        //- Construct null
        profilingPool(const IOobject&);


    //- Destructor
    ~profilingPool();

    
    // Member functions

        mapType& map()
        {
            return allInfo_;
        }

        const mapType& map() const
        {
            return allInfo_;
        }

        profilingStack& stack() 
        {
            return theStack_;
        }

        const profilingStack& stack() const
        {
            return theStack_;
        }

public:

    static void initprofiling(const IOobject&);

    static profilingInfo& getInfo(const string& name);

    static void remove(const profilingInfo& info);

    static void rememberTimer(const profilingInfo& info, clockTime& timer);

    virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
