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
    fineNumericFlux

Description
    Fine-level numeric flux class for density-based solvers

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    fineNumericFlux.H
    fineNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef fineNumericFlux_H
#define fineNumericFlux_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basicThermo.H"

#include "mgMeshLevel.H"
#include "mgFieldLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fineNumericFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Flux, class Limiter>
class fineNumericFlux
:
    public numericFluxBase<Flux>
{
    // Private data

        //- Reference to fine mesh level
        const mgMeshLevel& meshLevel_;

        //- Reference to fine fields level
        const mgFieldLevel& fieldLevel_;


        // Reference to primitive fields

            //- Static pressure
            const volScalarField& p_;

            //- Velocity
            const volVectorField& U_;

            //- Static temperature
            const volScalarField& T_;

            //- Reference to the thermophysicalModel
            basicThermo& thermo_;


        // Fluxes

            //- Density flux
            surfaceScalarField& rhoFlux_;

            //- Velocity flux
            surfaceVectorField& rhoUFlux_;

            //- Energy flux
            surfaceScalarField& rhoEFlux_;


        // Gradients: no longer stored.  HJ, 22/Apr/2016


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fineNumericFlux(const fineNumericFlux&);

        //- Disallow default bitwise assignment
        void operator=(const fineNumericFlux&);


        //- Return internal field of mass flux
        const scalarField& rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }

        //- Return access to internal field of mass flux
        scalarField& rhoFluxI()
        {
            return rhoFlux_.internalField();
        }

        //- Return internal field of momentum flux
        const vectorField& rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }

        //- Return access to internal field of momentum flux
        vectorField& rhoUFluxI()
        {
            return rhoUFlux_.internalField();
        }

         //- Return access to internal field of energy flux
        const scalarField& rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }

        //- Return access to internal field of energy flux
        scalarField& rhoEFluxI()
        {
            return rhoEFlux_.internalField();
        }


public:

    // Constructors

        //- Construct from components
        fineNumericFlux
        (
            const mgMeshLevel& meshLevel,
            const mgFieldLevel& fieldLevel,
            basicThermo& thermo
        );



    //- Destructor
    virtual ~fineNumericFlux()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return meshLevel_.mesh();
        }


        // Return fluxes

            //- Return density flux
            virtual const surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            //- Return velocity flux
            virtual const surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }

            //- Return energy flux
            virtual const surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }


       // Return residuals

            //- Return density equation residual
            virtual tmp<scalarField> rhoResidual() const
            {
                return fvc::div(rhoFlux_)().internalField();
            }

            //- Return momentum equation flux
            virtual tmp<vectorField> rhoUResidual() const
            {
                return fvc::div(rhoUFlux_)().internalField();
            }

            //- Return energy equation flux
            virtual tmp<scalarField> rhoEResidual() const
            {
                return fvc::div(rhoEFlux_)().internalField();
            }


        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "fineNumericFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
