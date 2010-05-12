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
    Foam::incompressible::RASModels::NonlinearKEShih

Description
    Shih's quadratic non-linear k-epsilon turbulence model for
    incompressible flows

SourceFiles
    NonlinearKEShih.C

\*---------------------------------------------------------------------------*/

#ifndef NonlinearKEShih_H
#define NonlinearKEShih_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class NonlinearKEShih Declaration
\*---------------------------------------------------------------------------*/

class NonlinearKEShih
:
    public RASModel
{
    // Private data

        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar alphak_;
        dimensionedScalar alphaEps_;
        dimensionedScalar A1_;
        dimensionedScalar A2_;
        dimensionedScalar Ctau1_;
        dimensionedScalar Ctau2_;
        dimensionedScalar Ctau3_;
        dimensionedScalar alphaKsi_;

        volScalarField k_;
        volScalarField epsilon_;

        volTensorField gradU_;
        volScalarField eta_;
        volScalarField ksi_;
        volScalarField Cmu_;
        volScalarField fEta_;

        volScalarField nut_;

        volSymmTensorField nonlinearStress_;

public:

    //- Runtime type information
    TypeName("NonlinearKEShih");

    // Constructors

        //- from components
        NonlinearKEShih
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~NonlinearKEShih()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphak_*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", alphaEps_*nut_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
