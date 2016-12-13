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
    quadraticDiff

Description
    Quadratic distance-based motion diffusion.

SourceFiles
    quadraticDiff.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticDiff_H
#define quadraticDiff_H

#include "linearDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class quadraticDiff Declaration
\*---------------------------------------------------------------------------*/

class quadraticDiff
:
    public linearDiff
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        quadraticDiff(const quadraticDiff&);

        //- Disallow default bitwise assignment
        void operator=(const quadraticDiff&);


public:

    //- Runtime type information
    TypeName("quadratic");


    // Constructors

        //- Construct from components
        quadraticDiff
        (
            const tetMotionSolver& mSolver
        );


    // Destructor

        virtual ~quadraticDiff();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
