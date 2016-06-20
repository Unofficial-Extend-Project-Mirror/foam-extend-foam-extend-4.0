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
    pisoFluid

Description
    icoFoam fluid solver

Author
    Zeljko Tukovic, FSB.  All rights reserved.

SourceFiles
    pisoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef pisoFluid_H
#define pisoFluid_H

#include "fluidSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidSolvers
{

/*---------------------------------------------------------------------------*\
                        Class pisoFluid Declaration
\*---------------------------------------------------------------------------*/

class pisoFluid
:
    public fluidSolver
{
    // Private data

        //- Velocity field
        volVectorField U_;

        //- Pressure field
        volScalarField p_;

        //- Pressure field
        volVectorField gradp_;

        //- Velocity gradient
        volTensorField gradU_;

        //- Flux field
        surfaceScalarField phi_;

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence_;

        //- Density
        dimensionedScalar rho_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pisoFluid(const pisoFluid&);

        //- Disallow default bitwise assignment
        void operator=(const pisoFluid&);


public:

    //- Runtime type information
    TypeName("pisoFluid");

    // Constructors

        //- Construct from components
        pisoFluid(const fvMesh& mesh);

//         //- Construct from components
//         pisoFluid
//         (
//             const word& type,
//             const fvMesh& mesh
//         );

    // Destructor

        virtual ~pisoFluid()
        {}


    // Member Functions

        // Access

            //- Return velocity field
            virtual const volVectorField& U() const;

            //- Return velocity field
            volVectorField& U()
            {
                return U_;
            }

            //- Return pressure field
            virtual const volScalarField& p() const;

            //- Return pressure field
            volScalarField& p()
            {
                return p_;
            }

            //- Return pressure gradient
            volVectorField& gradp()
            {
                return gradp_;
            }

            //- Return flux field
            surfaceScalarField& phi()
            {
                return phi_;
            }

//             //- Return kinematic viscosity
//             const dimensionedScalar& nu() const
//             {
//                 return nu_;
//             }

            //- Density
            const dimensionedScalar& rho()
            {
                return rho_;
            }

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID
            ) const;

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
