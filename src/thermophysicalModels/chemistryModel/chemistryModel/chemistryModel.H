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
    Foam::chemistryModel

Description
    Foam::chemistryModel

SourceFiles
    chemistryModelI.H
    chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryModel_H
#define chemistryModel_H

#include "hCombustionThermo.H"
#include "reactingMixture.H"
#include "Matrix.H"
#include "ODE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class chemistrySolver;

/*---------------------------------------------------------------------------*\
                           Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

class chemistryModel
:
    public ODE
{

public:

    typedef reactingMixture::reaction reaction;
    typedef reactingMixture::reactionThermo reactionThermo;


private:

    // Private data

        const hCombustionThermo& thermo_;

        PtrList<volScalarField>& Y_;
        const volScalarField& rho_;

        IOdictionary chemistryProperties_;

        //- Do the chemistry or not.
        Switch chemistry_;

        //- The reactions
        const PtrList<reaction>& reactions_;

        //- thermodynamic data of the species
        const PtrList<reactionThermo>& specieThermo_;

        //- The number of species
        label Ns_;

        //- The number of reactions
        label Nr_;

        //- ODE coefficients
        scalarField coeffs_;

        //- Chemistry solver
        autoPtr<chemistrySolver> solver_;

        //- the chemical source term
        PtrList<scalarField> RR_;

        //- Latest estimation of integration step
        scalarField deltaTChem_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const chemistryModel&);


public:

    // Constructors

        //- Construct from components
        chemistryModel
        (
            hCombustionThermo& thermo,
            const volScalarField& rho
        );


    // Destructor

        virtual ~chemistryModel();


    // Member Functions

        //- Return the thermodynamics model
        const hCombustionThermo& thermo() const
        {
            return thermo_;
        }

        //- Access the fractions of the species
        const PtrList<volScalarField>& Y() const
        {
            return Y_;
        }

        //- Access the density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the chemistry properties
        inline const IOdictionary& chemistryProperties() const
        {
            return chemistryProperties_;
        }

        //- Do the chemistry or not.
        inline Switch chemistry() const
        {
            return chemistry_;
        }

        inline const PtrList<reaction>& reactions() const
        {
            return reactions_;
        }

        inline const PtrList<reactionThermo>& specieThermo() const
        {
            return specieThermo_;
        }

        //- The number of species
        inline label Ns() const
        {
            return Ns_;
        }

        //- The number of reactions
        inline label Nr() const
        {
            return Nr_;
        }

        // Return the chemisty solver
        inline const chemistrySolver& solver() const
        {
            return solver_;
        }

        //- Return the chemical source terms
        inline tmp<volScalarField> RR(const label i) const;

        //- Return the latest estimation of integration step
        const scalarField& deltaTChem() const
        {
            return deltaTChem_;
        }

        //- The chemical time scale
        tmp<volScalarField> tc() const;

        //- dc/dt = omega, rate of change in concentration, for each species
        scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- return the reaction rate for reaction r and the reference
        //  species and charateristic times
        scalar omega
        (
            const reaction& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- The number of ODE to solve
        label nEqns() const;

        //- Return access to coefficients
        virtual scalarField& coeffs()
        {
            return coeffs_;
        }

        //- Return reference to coefficients
        virtual const scalarField& coeffs() const
        {
            return coeffs_;
        }

        void derivatives
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt
        ) const;

        void jacobian
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt,
            Matrix<scalar>& dfdc
        ) const;

        //- Update ODE after the solution, advancing by delta
        virtual void update(const scalar delta)
        {}

        //- Calculates the reaction rates
        void calculate();

        //- Solve the reaction system for the given start time and time-step
        // and return the characteristic time
        scalar solve
        (
            const scalar t0,
            const scalar deltaT
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
