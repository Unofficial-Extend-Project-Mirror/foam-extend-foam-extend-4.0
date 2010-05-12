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
    Foam::incompressible::RASModels::LRR

Description
    Launder, Reece and Rodi Reynolds-stress turbulence model for
    incompressible flows.

    The default model coefficients correspond to the following:
    @verbatim
        LRRCoeffs
        {
            Cmu         0.09;
            Clrr1       1.8;
            Clrr2       0.6;
            C1          1.44;
            C2          1.92;
            Cs          0.25;
            Ceps        0.15;
            alphah      1.0;    // only for compressible
            alphaEps    0.76923;
            alphaR      1.22;   // only for compressible
            couplingFactor  0.0;    // only for incompressible
        }
    @endverbatim

SourceFiles
    LRR.C

\*---------------------------------------------------------------------------*/

#ifndef LRR_H
#define LRR_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LRR Declaration
\*---------------------------------------------------------------------------*/

class LRR
:
    public RASModel
{
    // Private data

        dimensionedScalar Cmu_;

        dimensionedScalar Clrr1_;
        dimensionedScalar Clrr2_;

        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar Cs_;
        dimensionedScalar Ceps_;
        dimensionedScalar alphaEps_;

        dimensionedScalar couplingFactor_;

        volSymmTensorField R_;
        volScalarField k_;
        volScalarField epsilon_;
        volScalarField nut_;


public:

    //- Runtime type information
    TypeName("LRR");

    // Constructors

        //- from components
        LRR
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~LRR()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for R
        tmp<volScalarField> DREff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DREff", nut_ + nu())
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
        tmp<volSymmTensorField> R() const
        {
            return R_;
        }

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
