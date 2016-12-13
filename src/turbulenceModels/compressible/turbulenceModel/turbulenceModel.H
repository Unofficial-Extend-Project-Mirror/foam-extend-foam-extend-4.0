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

Namespace
    Foam::compressible::turbulenceModels

Description
    Namespace for compressible turbulence turbulence models.


Class
    Foam::compressible::turbulenceModel

Description
    Abstract base class for compressible turbulence models
    (RAS, LES and laminar).


SourceFiles
    turbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleturbulenceModel_H
#define compressibleturbulenceModel_H

#include "objectRegistry.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "basicThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace compressible
{

/*---------------------------------------------------------------------------*\
                       Class turbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModel
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volScalarField& rho_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        const basicThermo& thermophysicalModel_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        turbulenceModel(const turbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModel&);


public:

    //- Runtime type information
    TypeName("turbulenceModel");


    // Declare run-time constructor selection table

#ifndef SWIG
        declareRunTimeNewSelectionTable
        (
            autoPtr,
            turbulenceModel,
            turbulenceModel,
            (
                const volScalarField& rho,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const basicThermo& thermoPhysicalModel
            ),
            (rho, U, phi, thermoPhysicalModel)
        );
#endif


    // Constructors

        //- Construct from components
        turbulenceModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<compressible::turbulenceModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    //- Destructor
    virtual ~turbulenceModel()
    {}


    // Member Functions

        //- Access function to density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Access function to velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Access function to thermophysical model
        const basicThermo& thermo() const
        {
             return thermophysicalModel_;
        }

        //- Return the laminar viscosity
        const volScalarField& mu() const
        {
            return thermophysicalModel_.mu();
        }

        //- Return the laminar thermal conductivity
        const volScalarField& alpha() const
        {
            return thermophysicalModel_.alpha();
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> muEff() const = 0;

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
        virtual tmp<fvVectorMatrix> divDevRhoReff() const = 0;

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
