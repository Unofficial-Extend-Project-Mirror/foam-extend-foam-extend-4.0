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
    XPP_DE

Description
    eXtended Pom-Pom (XPP) double equation non linear viscoelastic fluid model
    (Verbeeten, W. M. H., Computational Polymer Melt Rheology, 2001).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    XPP_DE.C

\*---------------------------------------------------------------------------*/

#ifndef XPP_DE_H
#define XPP_DE_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class XPP_DE Declaration
\*---------------------------------------------------------------------------*/

class XPP_DE
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

            //- Material parameter, defining the amount of anisotropy
            dimensionedScalar alpha_;

            //- Relaxation time of the backbone tube orientation
            dimensionedScalar lambdaOb_;

            //- Relaxation time for the stretch
            dimensionedScalar lambdaOs_;

            //- Amount of arms at the end of a backbone
            dimensionedScalar q_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        XPP_DE(const XPP_DE&);

        //- Disallow default bitwise assignment
        void operator=(const XPP_DE&);


public:

    //- Runtime type information
    TypeName("XPP-DE");

    // Constructors

        //- Construct from components
        XPP_DE
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~XPP_DE()
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
