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
    Foam::localTimeStep

Description
    The time-step is computed locally so that the local Courant number
    does not exceed the specified value.

Author
    Oliver Borm
    Aleksandar Jemcov

SourceFiles
    localTimeStep.C

\*---------------------------------------------------------------------------*/

#ifndef localTimeStep_H
#define localTimeStep_H

#include "fvCFD.H"
#include "basicThermo.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class localTimeStep Declaration
\*---------------------------------------------------------------------------*/

class localTimeStep
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        // Reference to primitive fields

            //- Velocity field
            const volVectorField& U_;

           //- Reference to the thermophysicalModel
           const basicThermo& thermo_;

           //- Reference to the turbulenceModel
           const compressible::turbulenceModel& turbModel_;


       //- CFL number limited local time step field
       volScalarField CoDeltaT_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        localTimeStep(const localTimeStep&);

        //- Disallow default bitwise assignment
        void operator=(const localTimeStep&);


public:

    // Constructors

        //- Construct from mesh and Istream
        localTimeStep
        (
            const volVectorField& U,
            const basicThermo& thermo,
            compressible::turbulenceModel& turbModel
        );


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the Courant-number limited time-step
        const volScalarField& CoDeltaT() const
        {
            return CoDeltaT_;
        }

        //- Update the local deltaT field
        // and adjust it to the minimum value optionally
        void update(const scalar maxCo, const bool adjustTimeStep = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
