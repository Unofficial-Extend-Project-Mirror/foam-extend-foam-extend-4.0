// The FOAM Project // File: elasticPlastic.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   elasticPlastic
   \\  /           | Family: rheologyLaw
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    nonLinearElasticPlastic

DESCRIPTION
    Linear-elastic-plastic rheology

C++
    listOfSourceFiles
        nonLinearElasticPlastic.C
    endListOfSourceFiles

AUTHOR
    Hrvoje Jasak.
*/
// ------------------------------------------------------------------------- //

#ifndef nonLinearElasticPlastic_H
#define nonLinearElasticPlastic_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class nonLinearElasticPlastic
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

        //- Yield stress
        dimensionedScalar sigmaY_;

        //- Plastic modulus
        dimensionedScalar Ep_;

        //- Strentgh
        dimensionedScalar matStrength_;

        //- b coefficient of power law
        dimensionedScalar bCf_;

        //- n coefficient of power law
        dimensionedScalar nCf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        nonLinearElasticPlastic(const nonLinearElasticPlastic&);

        //- Disallow default bitwise assignment
        void operator=(const nonLinearElasticPlastic&);


public:

    //- Runtime type information
    TypeName("nonLinearElasticPlastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        nonLinearElasticPlastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~nonLinearElasticPlastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const;

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

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep() const;

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep(const volScalarField& sigmaEq) const;

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
