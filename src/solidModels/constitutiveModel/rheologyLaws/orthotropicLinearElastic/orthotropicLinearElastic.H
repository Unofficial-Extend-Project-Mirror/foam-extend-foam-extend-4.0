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
    orthotropicLinearElastic

Description
    Orthotropic linear elastic rheology
    with general material directions

    9 independent material parameters
    Ex,Ey,Ez,nuxy,nuyz,nuzx,Gxy,Gyz,Gzx

    Properties are specified for local coordinate system and are then
    rotated to the global cartesian coordinates based on the specified
    material direction fields.

    Please cite:
    Cardiff P, Karac A & Ivankovic A, A Large Strain Finite Volume Method for
    Orthotropic Bodies with General Material Orientations, Computer Methods
    in Applied Mechanics & Engineering, Sep 2013,
    http://dx.doi.org/10.1016/j.cma.2013.09.008

SourceFiles
    orthotropicLinearElastic.C

Author
    Philip Cardiff UCD
\*---------------------------------------------------------------------------*/

#ifndef orthotropicLinearElastic_H
#define orthotropicLinearElastic_H

#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class orthotropicLinearElastic Declaration
\*---------------------------------------------------------------------------*/

class orthotropicLinearElastic
:
    public rheologyLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Modulus of elasticity
        dimensionedScalar Ex_;
        dimensionedScalar Ey_;
        dimensionedScalar Ez_;
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nuxy_;
        dimensionedScalar nuyz_;
        dimensionedScalar nuzx_;
        dimensionedScalar nuyx_;
        dimensionedScalar nuzy_;
        dimensionedScalar nuxz_;
        dimensionedScalar nu_;

        //- shear modulus
        dimensionedScalar Gxy_;
        dimensionedScalar Gyz_;
        dimensionedScalar Gzx_;

        //- coefficients
        dimensionedScalar J_;
        dimensionedScalar A11_;
        dimensionedScalar A22_;
        dimensionedScalar A33_;
        dimensionedScalar A12_;
        dimensionedScalar A31_;
        dimensionedScalar A23_;
        dimensionedScalar A44_;
        dimensionedScalar A55_;
        dimensionedScalar A66_;

        //- stiffness fourth order tensor field
        //- I probably don't need to store C but
        //- then I have to make the rotation field every time
        //- C is asked for
        volSymmTensor4thOrderField C_;

        //- material directions tensor field
        //- (xx xy xz) is matDir1
        //- (yx yy yz) is matDir2
        //- (zx zy zz) is matDir3
        volTensorField matDir_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        orthotropicLinearElastic(const orthotropicLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const orthotropicLinearElastic&);

public:

    //- Runtime type information
    TypeName("orthotropicLinearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        orthotropicLinearElastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~orthotropicLinearElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

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

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const
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

        //- Return modulus of plasticity
        virtual tmp<volScalarField> Ep() const;

        virtual tmp<volScalarField> Ep(const volScalarField& sigmaEq) const
        {
            return Ep();
        }

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const;

        //- Return yield stress given cellID
        virtual scalar sigmaY
        (const scalar epsilonPEq, const label cellID) const;

        //- orthotropic linear elastic properties
        virtual tmp<volDiagTensorField> K() const;

        //- elastic constitutive fourth order tensor
        virtual tmp<volSymmTensor4thOrderField> C() const;

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
