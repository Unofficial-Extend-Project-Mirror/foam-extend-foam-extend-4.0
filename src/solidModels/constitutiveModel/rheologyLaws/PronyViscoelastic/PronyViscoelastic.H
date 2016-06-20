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

Class
    PronyViscoelastic

Description
    Linear elastic rheology

SourceFiles
    PronyViscoelastic.C

\*---------------------------------------------------------------------------*/

#ifndef PronyViscoelastic_H
#define PronyViscoelastic_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PronyViscoelastic Declaration
\*---------------------------------------------------------------------------*/

class PronyViscoelastic
:
    public rheologyLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Linear spring constants of generalised Maxwell model
        Field<scalar> k_;

        //- Dimension set for k field
        dimensionSet kDimensions_;

        //- Relaxation times of generalised Maxwell model
        Field<scalar> tau_;

        //- Dimension set for tau field
        dimensionSet tauDimensions_;

        //- Poisson's ratio
        dimensionedScalar nu_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PronyViscoelastic(const PronyViscoelastic&);

        //- Disallow default bitwise assignment
        void operator=(const PronyViscoelastic&);


public:

    //- Runtime type information
    TypeName("PronyViscoelastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        PronyViscoelastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~PronyViscoelastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const
        {
            return rho(0);
        }

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const
        {
            return E(0);
        }

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const
        {
            return nu(0);
        }

        //- Return density
        virtual tmp<volScalarField> rho(scalar t) const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(scalar t) const;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const;

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const;

        //- Return modulus of plasticity
        virtual tmp<volScalarField> Ep() const;

        virtual tmp<volScalarField> Ep(const volScalarField& sigmaEq) const
        {
            return Ep();
        }

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const;

        //- Return yield stress
        virtual scalar sigmaY(const scalar epsilonPEq, const label cellID) const
        {
            return GREAT;
        }

        //- If effects are active
        virtual bool viscoActive() const
        {
            return true;
        };

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
