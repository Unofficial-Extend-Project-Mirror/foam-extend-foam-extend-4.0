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
    Foam::compressible::LESModels

Description
    Namespace for compressible LES models.


Class
    Foam::compressible::LESModel

Description
    Class for all compressible flow LES SGS models.

    This class defines the basic interface for a compressible flow SGS model,
    and encapsulates data of value to all possible models. In particular
    this includes references to all the dependent fields (rho, U, phi),
    the physical viscosity mu, and the turbulenceProperties dictionary
    which contains the model selection and model coefficients.

SourceFiles
    LESModel.C
    newLESModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleLESModel_H
#define compressibleLESModel_H

#include "LESdelta.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class LESModel Declaration
\*---------------------------------------------------------------------------*/

class LESModel
:
    public IOdictionary
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;


private:

    // Private data

        const volScalarField& rho_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        const basicThermo& thermoPhysicalModel_;

        Switch printCoeffs_;
        dictionary coeffDict_;

        dimensionedScalar k0_;

        autoPtr<LESdelta> delta_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LESModel(const LESModel&);
        LESModel& operator=(const LESModel&);


protected:

    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


public:

    //- Runtime type information
    TypeName("LESModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESModel,
            dictionary,
            (
                const volScalarField& rho,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const basicThermo& thermoPhysicalModel
            ),
            (rho, U, phi, thermoPhysicalModel)
        );


    // Constructors

        //- Construct from components
        LESModel
        (
            const word& type,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<LESModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    // Destructor

        virtual ~LESModel()
        {}


    // Member Functions

        // Access

            //- Access function to the density field
            inline const volScalarField& rho() const
            {
                return rho_;
            }

            //- Access function to velocity field
            inline const volVectorField& U() const
            {
                return U_;
            }

            //- Access function to flux field
            inline const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Access function to the thermophysical properties model
            inline const basicThermo& thermo() const
            {
                return thermoPhysicalModel_;
            }

            //- Access the dictionary which provides info. about choice of
            //  models, and all related data (particularly model coefficients).
            inline dictionary& coeffDict()
            {
                return coeffDict_;
            }

            //- Access function to filter width
            inline const volScalarField& delta() const
            {
                return delta_();
            }

            //- Return the value of k0 which k is not allowed to be less than
            const dimensionedScalar& k0() const
            {
                return k0_;
            }

            //- Allow k0 to be changed
            dimensionedScalar& k0()
            {
                return k0_;
            }

            //- Access function to laminar viscosity
            tmp<volScalarField> mu() const
            {
                return thermoPhysicalModel_.mu();
            }

            //- Access function to laminar thermal conductivity
            tmp<volScalarField> alpha() const
            {
                return thermoPhysicalModel_.alpha();
            }


        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const = 0;

        //- Return the SGS turbulent dissipation.
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> muSgs() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("muEff", muSgs() + mu())
            );
        }

        //- Return the SGS thermal conductivity.
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const = 0;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoBeff() const = 0;

        //- Returns div(rho*dev(B)).
        // This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevRhoBeff(volVectorField& U) const = 0;

        //- Correct Eddy-Viscosity and related properties.
        //  This calls correct(const tmp<volTensorField>& gradU) by supplying
        //  gradU calculated locally.
        void correct();

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read turbulenceProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
