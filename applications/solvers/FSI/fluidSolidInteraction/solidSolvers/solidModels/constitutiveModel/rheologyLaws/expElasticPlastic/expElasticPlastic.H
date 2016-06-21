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
    expElasticPlastic

Description
    Non-linear plastic constitutive law.

SourceFiles
    expElasticPlastic.C

Author
    Zeljko Tukovic, FSB Zagreb

\*---------------------------------------------------------------------------*/

#ifndef expElasticPlastic_H
#define expElasticPlastic_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class expElasticPlastic
:
    public rheologyLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Modulus of elasticity
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Initial Yield stress
        dimensionedScalar sigmaY_;

        //- Plastic modulus
        dimensionedScalar Ep_;

        //- Saturation Yield stress
        dimensionedScalar sigmaYInf_;

        //- Saturation exponent
        dimensionedScalar delta_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        expElasticPlastic(const expElasticPlastic&);

        //- Disallow default bitwise assignment
        void operator=(const expElasticPlastic&);


public:

    //- Runtime type information
    TypeName("expElasticPlastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        expElasticPlastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~expElasticPlastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const;

        //- Return density
        virtual tmp<volScalarField> rho(scalar t) const
        {
            return rho();
        }

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(scalar t) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const
        {
            return nu();
        }

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const
        {
            notImplemented(type() + "::J(scalar t)");

            return 1.0/E(t);
        }

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const;

        //- Return yield stress
        virtual scalar sigmaY
        (
            const scalar epsilonPEq,
            const label cellID
        ) const;

        //- Return derivative of the yield stress curve
        virtual scalar dSigmaY
        (
            const scalar epsilonPEq,
            const label cellID
        ) const;

        //- If the rheologyLaw needs plasticity procedures
        virtual bool plasticityModelNeeded() const
        {
            return true;
        };

        //- If the plasticity law is nonlinear
        // This lets the plasticity model decide if an iterative method is
        // required
        virtual bool nonLinearPlasticity() const
        {
            return true;
        };

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep() const;

        virtual tmp<volScalarField> Ep(const volScalarField& sigmaEq) const
        {
            return Ep();
        }

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
