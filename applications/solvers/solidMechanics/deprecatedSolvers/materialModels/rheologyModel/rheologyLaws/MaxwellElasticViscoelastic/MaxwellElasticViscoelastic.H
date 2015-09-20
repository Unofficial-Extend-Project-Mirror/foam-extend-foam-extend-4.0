/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    MaxwellElasticViscoelastic

Description
    Linear elastic rheology

SourceFiles
    MaxwellElasticViscoelastic.C

\*---------------------------------------------------------------------------*/

#ifndef MaxwellElasticViscoelastic_H
#define MaxwellElasticViscoelastic_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class MaxwellElasticViscoelastic Declaration
\*---------------------------------------------------------------------------*/

class MaxwellElasticViscoelastic
:
    public rheologyLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Linear spring constant
        dimensionedScalar k_;

        //- Coefficient of viscosity
        dimensionedScalar eta_;

        //- Poisson's ratio
        dimensionedScalar nu_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MaxwellElasticViscoelastic(const MaxwellElasticViscoelastic&);

        //- Disallow default bitwise assignment
        void operator=(const MaxwellElasticViscoelastic&);


public:

    //- Runtime type information
    TypeName("MaxwellElasticViscoelastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        MaxwellElasticViscoelastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~MaxwellElasticViscoelastic();


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

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
