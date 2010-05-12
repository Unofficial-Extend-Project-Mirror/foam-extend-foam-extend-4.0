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
    Foam::incompressible::LESModels::SpalartAllmaras

Description
    SpalartAllmaras for incompressible flows

SourceFiles
    SpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmaras_H
#define SpalartAllmaras_H

#include "LESModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public LESModel
{
    // Private data

        dimensionedScalar alphaNut_;

        dimensionedScalar Cb1_;
        dimensionedScalar Cb2_;
        dimensionedScalar Cv1_;
        dimensionedScalar Cv2_;
        dimensionedScalar CDES_;
        dimensionedScalar ck_;
        dimensionedScalar kappa_;
        dimensionedScalar Cw1_;
        dimensionedScalar Cw2_;
        dimensionedScalar Cw3_;


    // Private member functions

        tmp<volScalarField> fv1() const;
        tmp<volScalarField> fv2() const;
        tmp<volScalarField> fv3() const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;

        // Disallow default bitwise copy construct and assignment
        SpalartAllmaras(const SpalartAllmaras&);
        SpalartAllmaras& operator=(const SpalartAllmaras&);

        volScalarField nuTilda_;
        volScalarField dTilda_;
        volScalarField nuSgs_;


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- Constructor from components
        SpalartAllmaras
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~SpalartAllmaras()
        {}


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return sqr(nuSgs()/ck_/dTilda_);
        }

        //- Return sub-grid disipation rate
        tmp<volScalarField> epsilon() const;

        tmp<volScalarField> nuTilda() const
        {
            return nuTilda_;
        }

        //- Return SGS viscosity
        tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        tmp<volSymmTensorField> devBeff() const;

        //- Return the deviatoric part of the divergence of Beff
        //  i.e. the additional term in the filtered NSE.
        tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Correct nuTilda and related properties
        void correct(const tmp<volTensorField>& gradU);

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
