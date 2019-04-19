/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    sixDOFBodies

Description
    6-DOF solver for multiple bodies

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

SourceFiles
    sixDOFBodies.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFBodies_H
#define sixDOFBodies_H

#include "foamTime.H"
#include "IOdictionary.H"
#include "sixDOFODE.H"
#include "ODESolver.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sixDOFBodies Declaration
\*---------------------------------------------------------------------------*/

class sixDOFBodies
:
    public IOdictionary
{
    // Private data

        //- Reference to time
        const Time& runTime_;

        // Pointer list of Foam::sixDOFODE objects
        PtrList<sixDOFODE> odes_;

        // Pointer list of Foam::sixDOFODE solvers
        PtrList<ODESolver> solvers_;

        // Name list of solved bodies
        const wordList names_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDOFBodies(const sixDOFBodies&);

        //- Disallow default bitwise assignment
        void operator=(const sixDOFBodies&);

        //- Set bodies
        void setBodies();


public:

    // Constructors

        //- Construct from Time
        sixDOFBodies(const Time& runTime);


    // Destructor

        virtual ~sixDOFBodies()
        {}


    // Member Functions

        //- Return list of names
        const wordList& names() const;

        //- Return list of bodies
        const PtrList<sixDOFODE>& operator()() const;

        //- Solve ODEs and update position
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
