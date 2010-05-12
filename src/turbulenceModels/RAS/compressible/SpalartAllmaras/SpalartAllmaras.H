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
    Foam::compressible::RASModels::SpalartAllmaras

Description
    Spalart-Allmaras one-eqn mixing-length model for compressible
    external flows.

    Reference:
    @verbatim
        "A One-Equation Turbulence Model for Aerodynamic Flows"
        P.R. Spalart,
        S.R. Allmaras,
        La Recherche A´rospatiale, No. 1, 1994, pp. 5–21.

        Extended according to:

        "An Unstructured Grid Generation and Adaptive Solution Technique
        for High Reynolds Number Compressible Flows"
        G.A. Ashford,
        Ph.D. thesis, University of Michigan, 1996.
    @endverbatim

    The default model coefficients correspond to the following:
    @verbatim
        SpalartAllmarasCoeffs
        {
            Cb1         0.1355;
            Cb2         0.622;
            Cw2         0.3;
            Cw3         2.0;
            Cv1         7.1;
            Cv2         5.0;
            alphaNut    1.5;
            alphah      1.0;    // only for compressible
        }
    @endverbatim

SourceFiles
    SpalartAllmaras.C
    SpalartAllmarasCorrect.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSpalartAllmaras_H
#define combressibleSpalartAllmaras_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public RASModel
{
    // Private data

        dimensionedScalar alphaNut_;
        dimensionedScalar alphah_;

        dimensionedScalar Cb1_;
        dimensionedScalar Cb2_;
        dimensionedScalar Cw1_;
        dimensionedScalar Cw2_;
        dimensionedScalar Cw3_;
        dimensionedScalar Cv1_;
        dimensionedScalar Cv2_;

        volScalarField nuTilda_;
        volScalarField mut_;

        wallDist d_;


    // Private member functions

        tmp<volScalarField> chi() const;
        tmp<volScalarField> fv1(const volScalarField& chi) const;
        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fv3
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- from components
        SpalartAllmaras
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            basicThermo& thermophysicalModel
        );


    // Destructor

        ~SpalartAllmaras()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DnuTildaEff",
                    alphaNut_*rho_*nuTilda_ + mu()
                )
            );
        }
	      //- Return the effective turbulent thermal diffusivity
        tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("alphaEff", alphah_*mut_ + alpha())
            );
        }

        //- Return the turbulence kinetic energy
        tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
                )
            );
        }

        //- Return the turbulence kinetic energy dissipation rate
        tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epslion",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
                )
            );
        }

        //- Return the Reynolds stress tensor
        tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
