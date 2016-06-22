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
    multiMaterial

Description
    Zoned multi-material rheology controlled by an indicator field.

SourceFiles
    multiMaterial.C

\*---------------------------------------------------------------------------*/

#ifndef multiMaterial_H
#define multiMaterial_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class multiMaterial Declaration
\*---------------------------------------------------------------------------*/

class multiMaterial
:
    public rheologyLaw,
    public PtrList<rheologyLaw>
{
    // Private data

        //- Material indicator field
        volScalarField materials_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterial(const multiMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterial&);


        //- Calculate indicator field given index
        tmp<scalarField> indicator(const label index) const;

        //- Calculate indicator field given index
        scalar indicator(const label index, const label cellID) const;

public:

    //- Runtime type information
    TypeName("multiMaterial");

    // Static data members


    // Constructors

        //- Construct from dictionary
        multiMaterial
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMaterial();


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
        virtual tmp<volScalarField> nu(scalar) const
        {
            return nu();
        }

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const
        {
            notImplemented(type() + "::J(scalar t)");

            return 1.0/E(t);
        }

        //- If the rheologyLaw needs plasticity procedures
        virtual bool plasticityModelNeeded() const;

        //- Return modulus of plasticity
        virtual tmp<volScalarField> Ep() const;

        virtual tmp<volScalarField> Ep(const volScalarField& epsEq) const;

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const;

        //- Return yield stress
        virtual scalar sigmaY
        (const scalar epsilonPEq, const label cellID) const;

        //- Orthotropic linear elastic properties
        virtual tmp<volDiagTensorField> K() const;
        virtual tmp<volSymmTensor4thOrderField> C() const;

        //- Correct the rheological model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
