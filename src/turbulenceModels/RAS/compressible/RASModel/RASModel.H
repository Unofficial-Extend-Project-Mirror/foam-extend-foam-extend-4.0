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

Namespace
    Foam::compressible::RASModels

Description
    Namespace for compressible RAS turbulence models.


Class
    Foam::compressible::RASModel

Description
    Abstract base class for turbulence models for compressible and combusting
    flows.

SourceFiles
    RASModel.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleRASModel_H
#define compressibleRASModel_H

#include "volFields.H"
#include "surfaceFields.H"
#include "nearWallDist.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class RASModel Declaration
\*---------------------------------------------------------------------------*/

class RASModel
:
    public IOdictionary
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volScalarField& rho_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        basicThermo& thermophysicalModel_;

        Switch turbulence_;
        Switch printCoeffs_;
        dictionary coeffDict_;

        dimensionedScalar kappa_;
        dimensionedScalar E_;

        scalar yPlusLam(const scalar kappa, const scalar E);
        scalar yPlusLam_;

        dimensionedScalar k0_;
        dimensionedScalar epsilon0_;
        dimensionedScalar epsilonSmall_;

        nearWallDist y_;


    // Protected member functions

        //- Print model coefficients
        virtual void printCoeffs();

        //- Return the laminar thermal conductivity
        const volScalarField& alpha() const
        {
            return thermophysicalModel_.alpha();
        }


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RASModel(const RASModel&);

        //- Disallow default bitwise assignment
        void operator=(const RASModel&);


public:

    //- Runtime type information
    TypeName("RASModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASModel,
            dictionary,
            (
                const volScalarField& rho,
                const volVectorField& U,
                const surfaceScalarField& phi,
                basicThermo& thermoPhysicalModel
            ),
            (rho, U, phi, thermoPhysicalModel)
        );


    // Constructors

        //- Construct from components
        RASModel
        (
            const word& type,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            basicThermo& thermoPhysicalModel
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<RASModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            basicThermo& thermoPhysicalModel
        );


    // Destructor

        virtual ~RASModel()
        {}


    // Member Functions

        // Access

            //- Return the value of k0 which k is not allowed to be less than
            const dimensionedScalar& k0() const
            {
                return k0_;
            }

            //- Return the value of epsilon0 which epsilon is not allowed to be
            //  less than
            const dimensionedScalar& epsilon0() const
            {
                return epsilon0_;
            }

            //- Return the value of epsilonSmall which is added to epsilon when
            //  calculating nut
            const dimensionedScalar& epsilonSmall() const
            {
                return epsilonSmall_;
            }


            //- Allow k0 to be changed
            dimensionedScalar& k0()
            {
                return k0_;
            }

            //- Allow epsilon0 to be changed
            dimensionedScalar& epsilon0()
            {
                return epsilon0_;
            }

            //- Allow epsilonSmall to be changed
            dimensionedScalar& epsilonSmall()
            {
                return epsilonSmall_;
            }


            //- Return kappa for use in wall-functions
            dimensionedScalar kappa() const
            {
                return kappa_;
            }

            //- Return E for use in wall-functions
            dimensionedScalar E() const
            {
                return E_;
            }

            //- Return y+ at the edge of the laminar sublayer
            //  for use in wall-functions
            scalar yPlusLam() const
            {
                return yPlusLam_;
            }

            //- Const access to the coefficients dictionary
            const dictionary& coeffDict() const
            {
                return coeffDict_;
            }


        //- Return the laminar viscosity
        const volScalarField& mu() const
        {
            return thermophysicalModel_.mu();
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("muEff", mut() + mu())
            );
        }

        //- Return the effective turbulent thermal diffusivity
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const = 0;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const = 0;

        //- Return yPlus for the given patch
        virtual tmp<scalarField> yPlus(const label patchI) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Read turbulenceProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
