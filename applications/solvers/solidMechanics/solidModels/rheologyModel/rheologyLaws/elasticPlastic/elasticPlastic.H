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
    elasticPlastic

DESCRIPTION
    Linear-elastic-plastic rheology

C++
    listOfSourceFiles
        elasticPlastic.C
    endListOfSourceFiles

AUTHOR
    Hrvoje Jasak.
*/
// ------------------------------------------------------------------------- //

#ifndef elasticPlastic_H
#define elasticPlastic_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class elasticPlastic
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


    // Private Member Functions

        //- Disallow default bitwise copy construct
        elasticPlastic(const elasticPlastic&);

        //- Disallow default bitwise assignment
        void operator=(const elasticPlastic&);


public:

    //- Runtime type information
    TypeName("elasticPlastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        elasticPlastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~elasticPlastic();


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
