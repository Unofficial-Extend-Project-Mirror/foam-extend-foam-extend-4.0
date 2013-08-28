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
    Foam::profilingStack

Description
    Stack of the profilingInfo-items that are currently used

SourceFiles
    profilingStack.C

\*---------------------------------------------------------------------------*/

#ifndef profilingStack_H
#define profilingStack_H

#include "LIFOStack.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
    class profilingInfo;
    class profilingPool;
    class Ostream;

/*---------------------------------------------------------------------------*\
                           Class profilingStack Declaration
\*---------------------------------------------------------------------------*/

class profilingStack
:
    private LIFOStack<profilingInfo *>
{

    // Private Member Functions

    //- Disallow default bitwise copy construct
    profilingStack(const profilingStack&);
    
    //- Disallow default bitwise assignment
    void operator=(const profilingStack&);

    //- remember the timers for the correct stack-output
    HashTable<clockTime *,label> timers_;
protected:
    void writeStackContents(Ostream &) const;

    void addTimer(const profilingInfo &info,clockTime &timer);
public:

    // Constructors

    //- Construct null
    profilingStack();

    // Destructor
    
    ~profilingStack();


    // Members that encapsulate the original stack-class
    profilingInfo &top() const;

    profilingInfo &bottom() const;
    
    bool empty() const;

    void push(profilingInfo &);

    profilingInfo &pop();

    friend class profilingPool;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
