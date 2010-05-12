/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    sixDOFbodies

Description
    6-DOF solver for multiple bodies

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

SourceFiles
    sixDOFbodies.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFbodies_H
#define sixDOFbodies_H

#include "Time.H"
#include "IOdictionary.H"
#include "sixDOFqODE.H"
#include "ODESolver.H"
#include "finiteRotation.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sixDOFbodies Declaration
\*---------------------------------------------------------------------------*/

class sixDOFbodies
:
    public IOdictionary
{
    // Private data

        //- Reference to time
        const Time& runTime_;

        // Pointer list of Foam::sixDOFqODE objects
        PtrList<sixDOFqODE> odes_;

        // Pointer list of Foam::sixDOFqODE solvers
        PtrList<ODESolver> solvers_;

        // Name list of solved bodies
        const wordList names_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDOFbodies(const sixDOFbodies&);

        //- Disallow default bitwise assignment
        void operator=(const sixDOFbodies&);

        //- Set bodies
        void setBodies();


public:

    // Constructors

        //- Construct from Time
        sixDOFbodies(const Time& runTime);


    // Destructor

        ~sixDOFbodies()
        {}


    // Member Functions

        //- Return list of names
        const wordList& names() const;

        //- Return list of bodies
        const PtrList<sixDOFqODE>& operator()() const;

        //- Solve ODEs and update position
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
