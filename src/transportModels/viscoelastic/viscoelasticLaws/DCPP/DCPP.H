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
    DCPP

Description
    Double Convected Pom-Pom (DCPP) viscoelastic fluid model
    (Nicolas Clemeur, Rulande P.G. Rutgers, Benoît Debbaut,
    Rheol Acta (2003) 42: 217-231).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    DCPP.C

\*---------------------------------------------------------------------------*/

#ifndef DCPP_H
#define DCPP_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DCPP Declaration
\*---------------------------------------------------------------------------*/

class DCPP
:
    public viscoelasticLaw
{
    // Private data

        //- Evolution of orientation
        volSymmTensorField S_;

        //- Evolution of the backbone stretch
        volScalarField Lambda_;

       //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- Identity tensor
        dimensionedSymmTensor I_;


        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- Material non linear parameter (0 ~ 0.2)
            dimensionedScalar zeta_;

            //- Relaxation time of the backbone tube orientation
            dimensionedScalar lambdaOb_;

            //- Relaxation time for the stretch
            dimensionedScalar lambdaOs_;

            //- Amount of arms at the end of a backbone
            dimensionedScalar q_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        DCPP(const DCPP&);

        //- Disallow default bitwise assignment
        void operator=(const DCPP&);


public:

    //- Runtime type information
    TypeName("DCPP");

    // Constructors

        //- Construct from components
        DCPP
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~DCPP()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

        //- Correct the viscoelastic stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
