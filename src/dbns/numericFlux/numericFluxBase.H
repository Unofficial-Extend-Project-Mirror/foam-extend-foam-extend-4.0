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
