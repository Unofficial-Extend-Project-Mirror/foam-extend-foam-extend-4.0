/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::incompressible::RASModels::LienCubicKE

Description
    Lien cubic non-linear low-Reynolds k-epsilon turbulence models for
    incompressible flows.

    This turbulence model is described in:
    \verbatim
        Lien, F.S., Chen, W.L. & Leschziner, M.A. (1996).
        Low-Reynolds-number eddy-viscosity modeling based on non-linear
        stress-strain/vorticity relations.
        Engineering Turbulence Modelling and Experiments 3, 91-100.
    \endverbatim

    Implemented according to the specification in:
    <a href=
    "http://personalpages.manchester.ac.uk/staff/david.d.apsley/specturb.pdf"
    >Apsley: Turbulence Models 2002</a>

    In addition to the low-Reynolds number damping functions support for
    wall-functions is also included to allow for low- and high-Reynolds number
    operation.

See Also
    Foam::incompressible::RASModels::ShihQuadraticKE

SourceFiles
    LienCubicKE.C

\*---------------------------------------------------------------------------*/

#ifndef LienCubicKE_H
#define LienCubicKE_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LienCubicKE Declaration
\*---------------------------------------------------------------------------*/

class LienCubicKE
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar A1_;
            dimensionedScalar A2_;
            dimensionedScalar Ctau1_;
            dimensionedScalar Ctau2_;
            dimensionedScalar Ctau3_;
            dimensionedScalar alphaKsi_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;

            volScalarField eta_;
            volScalarField ksi_;
            volScalarField Cmu_;
            volScalarField fEta_;
            volScalarField C5viscosity_;

            volScalarField nut_;

            volSymmTensorField nonlinearStress_;


public:

    //- Runtime type information
    TypeName("LienCubicKE");

    // Constructors

        //- Construct from components
        LienCubicKE
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~LienCubicKE()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
