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
    exponentialDiff

Description
    Exponential distance-based motion diffusion.

SourceFiles
    exponentialDiff.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialDiff_H
#define exponentialDiff_H

#include "linearDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class exponentialDiff Declaration
\*---------------------------------------------------------------------------*/

class exponentialDiff
:
    public linearDiff
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        exponentialDiff(const exponentialDiff&);

        //- Disallow default bitwise assignment
        void operator=(const exponentialDiff&);


public:

    //- Runtime type information
    TypeName("exponential");


    // Constructors

        //- Construct from components
        exponentialDiff
        (
            const tetDecompositionMotionSolver& mSolver
        );


    // Destructor

        virtual ~exponentialDiff();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
