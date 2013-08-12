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
    Foam::ProfilingStack

Description
    Stack of the ProfilingInfo-items that are currently used

SourceFiles
    ProfilingStack.C

\*---------------------------------------------------------------------------*/

#ifndef ProfilingStack_H
#define ProfilingStack_H

#include "LIFOStack.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
    class ProfilingInfo;
    class ProfilingPool;
    class Ostream;

/*---------------------------------------------------------------------------*\
                           Class ProfilingStack Declaration
\*---------------------------------------------------------------------------*/

class ProfilingStack
:
    private LIFOStack<ProfilingInfo *>
{

    // Private Member Functions

    //- Disallow default bitwise copy construct
    ProfilingStack(const ProfilingStack&);
    
    //- Disallow default bitwise assignment
    void operator=(const ProfilingStack&);

    //- remember the timers for the correct stack-output
    HashTable<clockTime *,label> timers_;
protected:
    void writeStackContents(Ostream &) const;

    void addTimer(const ProfilingInfo &info,clockTime &timer);
public:

    // Constructors

    //- Construct null
    ProfilingStack();

    // Destructor
    
    ~ProfilingStack();


    // Members that encapsulate the original stack-class
    ProfilingInfo &top() const;

    ProfilingInfo &bottom() const;
    
    bool empty() const;

    void push(ProfilingInfo &);

    ProfilingInfo &pop();

    friend class ProfilingPool;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
