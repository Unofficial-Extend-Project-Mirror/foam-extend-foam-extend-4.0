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
    FASNumericFlux

Description
    Full Approximation Storage Multigrid numeric flux

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    FASNumericFlux.H
    FASNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef FASNumericFlux_H
#define FASNumericFlux_H

#include "numericFluxBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class FASNumericFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Flux>
class FASNumericFlux
:
    public numericFluxBase<Flux>
{
    // Private data

        //- Reference to meshLevel
        const mgMeshLevel& meshLevel_;

        //- Reference to fieldLevel
        const mgFieldLevel& fieldLevel_;

        // Primitive fields

            //- Static pressure
            const scalarField p_;

            //- Velocity
            const vectorField U_;

            //- static Temperature
            const scalarField T_;


        // References to fine fluxes

            //- Density flux
            surfaceScalarField& fineRhoFlux_;

            //- momentum flux
            surfaceVectorField& fineRhoUFlux_;

            //- Energy flux
            surfaceScalarField& fineRhoEFlux_;


        // Coarse fluxes

            //- Density flux
            scalarField rhoFlux_;

            //- Momentum flux
            vectorField rhoUFlux_;

            //- Energy flux
            scalarField rhoEFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        FASNumericFlux(const FASNumericFlux&);

        //- Disallow default bitwise assignment
        void operator=(const FASNumericFlux&);


public:

    // Constructors

        //- Construct from components
        FASNumericFlux
        (
            const mgMeshLevel& meshLevel,
            const mgFieldLevel& fieldLevel
        );

    // Destructor

       virtual ~FASNumericFlux() {}


    // Member Functions

        //- Return reference to the mesh
        const fvMesh& mesh() const
        {
            return meshLevel_.mesh();
        }


        // Return fluxes

            //- Return density equation flux
            virtual const surfaceScalarField& rhoFlux() const
            {
                return fineRhoFlux_;
            }

            //- Return momentum equation flux
            virtual const surfaceVectorField& rhoUFlux() const
            {
                return fineRhoUFlux_;
            }

            //- Return energy equation flux
            virtual const surfaceScalarField& rhoEFlux() const
            {
                return fineRhoEFlux_;
            }


       // Return residuals

            //- Return density residual
            virtual tmp<scalarField> rhoResidual() const;

            //- Return momentum residual
            virtual tmp<vectorField> rhoUResidual() const;

            //- Return energy residual
            virtual tmp<scalarField> rhoEResidual() const;


        //- Update fluxes based on current state
        virtual void computeFlux();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "FASNumericFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
