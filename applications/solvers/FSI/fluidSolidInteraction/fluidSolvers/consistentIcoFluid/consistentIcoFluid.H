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
    consistentIcoFluid

Description
    icoFoam fluid solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    consistentIcoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef consistentIcoFluid_H
#define consistentIcoFluid_H

#include "icoFluid.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidSolvers
{

/*---------------------------------------------------------------------------*\
                        Class consistentIcoFluid Declaration
\*---------------------------------------------------------------------------*/

class consistentIcoFluid
:
    public icoFluid
{
    // Private data

        //- Face area vectors
        mutable surfaceVectorField* SfPtr_;

    // Private Member Functions

        //- Make face area vectors
        void makeSf() const;

        //- Make face area vectors
        void updateSf();

        //- Disallow default bitwise copy construct
        consistentIcoFluid(const consistentIcoFluid&);

        //- Disallow default bitwise assignment
        void operator=(const consistentIcoFluid&);


public:

    //- Runtime type information
    TypeName("consistentIcoFluid");

    // Constructors

        //- Construct from components
        consistentIcoFluid(const fvMesh& mesh);


    // Destructor

        virtual ~consistentIcoFluid()
        {}


    // Member Functions

        // Access

            //- Return face area vectors
            const surfaceVectorField& Sf() const;

            //- Return face area vectors
            surfaceVectorField& Sf();

        // Edit

            //- Evolve the fluid solver
            virtual void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
