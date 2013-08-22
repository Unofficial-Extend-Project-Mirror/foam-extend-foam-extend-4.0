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
    Foam::ProfilingInfo

Description
    Information needed for profiling

SourceFiles
    ProfilingInfo.C

\*---------------------------------------------------------------------------*/

#ifndef ProfilingInfo_H
#define ProfilingInfo_H

#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
    // class Istream;
    class Ostream;
    class ProfilingStack;
    class ProfilingPool;
    class ProfilingInfo;

    Ostream& operator<<(Ostream&, const ProfilingInfo&);

/*---------------------------------------------------------------------------*\
                           Class ProfilingInfo Declaration
\*---------------------------------------------------------------------------*/

class ProfilingInfo
{
    // Private data

    // nr of times this was called
    label calls_;

    // total time spent
    scalar totalTime_;

    // time spent in children
    scalar childTime_;

    // unique id to identify it
    label id_;

    // pointer to the parent object (if there is any)
    ProfilingInfo &parent_;

    // what this does
    string description_;

    // is this currently on the stack?
    bool onStack_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    ProfilingInfo(const ProfilingInfo&);
    
    //- Disallow default bitwise assignment
    void operator=(const ProfilingInfo&);

    // Static data members
    
    //- Counter that is used to generate the ids
    static label nextId_;

    //- get a new ID and update the counter
    static label getID();

    //- raise the next possible ID (to avoid ID-clashes during reading)
    void raiseID(label maxVal);

protected:
    
    void addedToStack() 
        { onStack_=true; }

    void removedFromStack() 
        { onStack_=false; }

    //- Construct null - only the master-element
    ProfilingInfo();

    void writeWithOffset(Ostream &os,bool offset=false,scalar time=0,scalar childTime=0) const; 

public:

    // Constructors
    
    //- Construct from components
    ProfilingInfo(ProfilingInfo &parent,const string &descr);

//     //- Construct from Istream
//     ProfilingInfo(Istream&);

    // Destructor

    ~ProfilingInfo();


    // Member Functions

    // Access
    
    label id() const
        { return id_; }

    label calls() const
        { return calls_; }

    scalar totalTime() const
        { return totalTime_; }

    scalar childTime() const
        { return childTime_; }

    bool onStack() const
        { return onStack_; }

    const string &description() const
        { return description_; }

    const ProfilingInfo &parent() const
        { return parent_; }

    //- Update it with a new timing information
    void update(scalar elapsedTime);

    friend class ProfilingStack;
    friend class ProfilingPool;

    // IOstream Operators

    //        friend Istream& operator>>(Istream&, ProfilingInfo&);
    friend Ostream& operator<<(Ostream&, const ProfilingInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
