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

Class
    numericFluxBase

Description
    Base numeric flux class for density-based solvers

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    numericFluxBase.H
    numericFluxBase.C

\*---------------------------------------------------------------------------*/

#ifndef numericFluxBase_H
#define numericFluxBase_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class numericFluxBase Declaration
\*---------------------------------------------------------------------------*/

template<class Flux>
class numericFluxBase
:
    public Flux
{
public:

    // Constructors

        //- Construct null
        numericFluxBase()
        {}


    //- Destructor
    virtual ~numericFluxBase()
    {}


    // Memeber functions

        // Return fluxes

            //- Return continuity equation flux
            virtual const surfaceScalarField& rhoFlux() const = 0;

            //- Return momentum equation flux
            virtual const surfaceVectorField& rhoUFlux() const = 0;

            //- Return energy equation flux
            virtual const surfaceScalarField& rhoEFlux() const = 0;


        // Return residuals

            //- Return mass residual
            virtual tmp<scalarField> rhoResidual() const = 0;

            //- Return momentum residual
            virtual tmp<vectorField> rhoUResidual() const = 0;

            //- Return energy residual
            virtual tmp<scalarField> rhoEResidual() const = 0;


        //- Update fluxes based on current state
        virtual void computeFlux() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
