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
    WhiteMetznerCarreauYasuda

Description
    White-Metzner Carreau-Yasuda Law non linear viscoelastic fluid model
    (Larson, R. G., Constitutive equations for polymer melts
    and solutions, 1988).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    WhiteMetznerCarreauYasuda.C

\*---------------------------------------------------------------------------*/

#ifndef WhiteMetznerCarreauYasuda_H
#define WhiteMetznerCarreauYasuda_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class WhiteMetznerCarreauYasuda Declaration
\*---------------------------------------------------------------------------*/

class WhiteMetznerCarreauYasuda
:
    public viscoelasticLaw
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;


        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- Relaxation time
            dimensionedScalar lambda_;

            //- Fitted parameter for etaP
            dimensionedScalar m_;

            //- Fitted parameter for lambda
            dimensionedScalar n_;

            //- Fitted parameter for etaP
            dimensionedScalar K_;

            //- Fitted parameter for lambda
            dimensionedScalar L_;

            //- Fitted parameter for etaP
            dimensionedScalar a_;

            //- Fitted parameter for lambda
            dimensionedScalar b_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WhiteMetznerCarreauYasuda(const WhiteMetznerCarreauYasuda&);

        //- Disallow default bitwise assignment
        void operator=(const WhiteMetznerCarreauYasuda&);


public:

    //- Runtime type information
    TypeName("WhiteMetznerCarreauYasuda");

    // Constructors

        //- Construct from components
        WhiteMetznerCarreauYasuda
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~WhiteMetznerCarreauYasuda()
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
