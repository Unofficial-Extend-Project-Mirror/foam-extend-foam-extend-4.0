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
    Foam::compressible::RASModels::realizableKE

Description
    Realizable k-epsilon turbulence model for compressible flows.

    Model described in the paper:
    @verbatim
        "A New k-epsilon Eddy Viscosity Model for High Reynolds Number
        Turbulent Flows"

        Tsan-Hsing Shih, William W. Liou, Aamir Shabbir, Zhigang Tang and
        Jiang Zhu

        Computers and Fluids Vol. 24, No. 3, pp. 227-238, 1995
    @endverbatim

    The default model coefficients correspond to the following:
    @verbatim
        realizableKE
        {
            Cmu         0.09;
            A0          4.0;
            C2          1.9;
            alphak      1.0;
            alphaEps    0.833333;
            alphah      1.0;    // only for compressible
        }
    @endverbatim

SourceFiles
    realizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef realizableKE_H
#define realizableKE_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class realizableKE Declaration
\*---------------------------------------------------------------------------*/

class realizableKE
:
    public RASModel
{
    // Private data

        dimensionedScalar Cmu_;
        dimensionedScalar A0_;
        dimensionedScalar C2_;
        dimensionedScalar alphak_;
        dimensionedScalar alphaEps_;
        dimensionedScalar alphah_;

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField mut_;

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU
        );

public:

    //- Runtime type information
    TypeName("realizableKE");

    // Constructors

        //- from components
        realizableKE
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            basicThermo& thermophysicalModel
        );


    // Destructor

        ~realizableKE(){}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphak_*mut_ + mu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", alphaEps_*mut_ + mu())
            );
        }

        //- Return the effective turbulent thermal diffusivity
        tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("alphaEff", alphah_*mut_ + alpha())
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
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
