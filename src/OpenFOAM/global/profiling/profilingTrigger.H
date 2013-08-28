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
    Foam::profilingTrigger

Description
    The object that does the actual measuring

SourceFiles
    profilingTrigger.C

\*---------------------------------------------------------------------------*/

#ifndef profilingTrigger_H
#define profilingTrigger_H

#include "clockTime.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class profilingInfo;
    class profilingPool;

/*---------------------------------------------------------------------------*\
                           Class profilingTrigger Declaration
\*---------------------------------------------------------------------------*/

class profilingTrigger
{
    // Private data

    clockTime clock_;

    profilingInfo &info_;

    bool running_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    profilingTrigger(const profilingTrigger&);

    //- Disallow default bitwise assignment
    void operator=(const profilingTrigger&);

protected:

    clockTime &clock() 
        { return clock_; }

    const profilingInfo &info() const
        { return info_; }

public:

    // Constructors

    profilingTrigger(const string &name);
        
    ~profilingTrigger();

    void stop();

    friend class profilingPool;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
