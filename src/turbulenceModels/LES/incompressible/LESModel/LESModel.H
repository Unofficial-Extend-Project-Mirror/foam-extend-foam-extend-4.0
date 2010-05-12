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
    Foam::incompressible::LESModels

Description
    Namespace for incompressible LES models.


Class
    Foam::incompressible::LESModel

Description
    Base class for all incompressible flow LES SGS models.

    This class defines the basic interface for an incompressible flow SGS
    model, and encapsulates data of value to all possible models. In
    particular this includes references to all the dependent fields (U,
    phi), the physical viscosity nu, and the turbulenceProperties
    dictionary which contains the model selection and model coefficients.

SourceFiles
    LESModel.C
    newIsoLESModel.C

\*---------------------------------------------------------------------------*/

#ifndef LESModel_H
#define LESModel_H

#include "LESdelta.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "incompressible/transportModel/transportModel.H"
#include "wallFvPatch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
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

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        transportModel& transport_;

        Switch printCoeffs_;
        dictionary coeffDict_;

        dimensionedScalar k0_;

        autoPtr<LESdelta> delta_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LESModel(const LESModel&);
        LESModel& operator=(const LESModel&);


protected:

    // Protected member functions

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
                const volVectorField& U,
                const surfaceScalarField& phi,
                transportModel& transport
            ),
            (U, phi, transport)
        );


    // Constructors

        //- Construct from components
        LESModel
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<LESModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        virtual ~LESModel()
        {}


    // Member Functions

        // Access

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

            //- Access function to incompressible transport model
            inline transportModel& transport() const
            {
                return transport_;
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

            //- Access function to incompressible viscosity
            inline const volScalarField& nu() const
            {
                return transport_.nu();
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


        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const = 0;

        //- Return the SGS turbulent dissipation.
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("nuEff", nuSgs() + nu())
            );
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const = 0;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const = 0;

        //- Returns div(dev(Beff)).
        //  This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const = 0;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Correct Eddy-Viscosity and related properties.
        //  This calls correct(const tmp<volTensorField>& gradU) by supplying
        //  gradU calculated locally.
        void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
