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
    Foam::incompressible::LESModels::laminar

Description
    LES model for laminar incompressible flow.
    It simply returns laminar properties.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef laminar_H
#define laminar_H

#include "LESModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public LESModel
{
    // Private member functions

        // Disallow default bitwise copy construct and assignment
        laminar(const laminar&);
        laminar& operator=(const laminar&);


public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Constructor from components
        laminar
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~laminar()
        {}


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const;

        //- Return sub-grid disipation rate
        tmp<volScalarField> epsilon() const;

        //- Return SGS viscosity
        tmp<volScalarField> nuSgs() const;

        //- Return the effective viscosity
        tmp<volScalarField> nuEff() const;

        //- Return the sub-grid stress tensor B.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const;

        //- Return the deviatoric part of the divergence of Beff
        //  i.e. the additional term in the filtered NSE.
        tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
