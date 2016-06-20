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
    uniformDiff

Description
    Uniform motion diffusion.

SourceFiles
    uniformDiff.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDiff_H
#define uniformDiff_H

#include "motionDiff.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class uniformDiff Declaration
\*---------------------------------------------------------------------------*/

class uniformDiff
:
    public motionDiff
{
private:

    // Private data

        elementScalarField motionGamma_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniformDiff(const uniformDiff&);

        //- Disallow default bitwise assignment
        void operator=(const uniformDiff&);


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from components
        uniformDiff
        (
            const tetMotionSolver& mSolver
        );


    // Destructor

        virtual ~uniformDiff();


    // Member Functions

        //- Return diffusivity field
        virtual const elementScalarField& motionGamma() const
        {
            return motionGamma_;
        }

        //- Correct the motion diffusivity
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
