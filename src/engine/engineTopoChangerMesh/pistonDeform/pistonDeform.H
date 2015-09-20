/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#ifndef pistonDeform_H
#define pistonDeform_H

#include "engineTopoChangerMesh.H"
#include "simpleEnginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class pistonDeform Declaration
\*---------------------------------------------------------------------------*/

class pistonDeform
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        simpleEnginePiston piston_;

        //- Piston Position
        scalar pistonPosition_;

        //- deckHeight
        scalar deckHeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pistonDeform(const pistonDeform&);

        //- Disallow default bitwise assignment
        void operator=(const pistonDeform&);

        //- Check if all patches exist, then calculate virtualPistonPosition,
        //  pistonPosition and deckHeight for the first time
        void checkAndCalculate();

        //- Add valve and piston zones and modifiers
        void addZonesAndModifiers()
        {}

public:

    //- Runtime type information
    TypeName("pistonDeform");


    // Constructors

        //- Construct from database
        explicit pistonDeform(const IOobject& io);


    // Destructor

        virtual ~pistonDeform();

    // Member Functions


        const scalar& pistonPosition() const
        {
            return pistonPosition_;
        }

        scalar& pistonPosition()
        {
            return pistonPosition_;
        }

        const scalar& deckHeight() const
        {
            return deckHeight_;
        }

        scalar& deckHeight()
        {
            return deckHeight_;
        }

        const simpleEnginePiston& piston() const
        {
            return piston_;
        }


        //- Move and morph
        virtual bool update();

        //- Set boundary velocities
        void setBoundaryVelocity(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
