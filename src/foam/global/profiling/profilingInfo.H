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

Class
    Foam::profilingInfo

Description
    Information needed for profiling

SourceFiles
    profilingInfo.C

\*---------------------------------------------------------------------------*/

#ifndef profilingInfo_H
#define profilingInfo_H

#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

// class Istream;
class Ostream;
class profilingStack;
class profilingPool;
class profilingInfo;

Ostream& operator<<(Ostream&, const profilingInfo&);

/*---------------------------------------------------------------------------*\
                           Class profilingInfo Declaration
\*---------------------------------------------------------------------------*/

class profilingInfo
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
    profilingInfo &parent_;

    // what this does
    string description_;

    // is this currently on the stack?
    bool onStack_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    profilingInfo(const profilingInfo&);

    //- Disallow default bitwise assignment
    void operator=(const profilingInfo&);

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
    profilingInfo();

    void writeWithOffset(Ostream &os,bool offset=false,scalar time=0,scalar childTime=0) const;

public:

    // Constructors

    //- Construct from components
    profilingInfo(profilingInfo &parent,const string &descr);

//     //- Construct from Istream
//     profilingInfo(Istream&);

    // Destructor

    ~profilingInfo();


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

    const profilingInfo &parent() const
        { return parent_; }

    //- Update it with a new timing information
    void update(scalar elapsedTime);

    friend class profilingStack;
    friend class profilingPool;

    // IOstream Operators

    //        friend Istream& operator>>(Istream&, profilingInfo&);
    friend Ostream& operator<<(Ostream&, const profilingInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
