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
    Feta_PTT

Description
    Fixed eta Phan-Thin-Tanner non linear viscoelastic fluid model
    (Verbeeten, W. M. H., Computational Polymer Melt Rheology, 2001).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    Feta_PTT.C

\*---------------------------------------------------------------------------*/

#ifndef Feta_PTT_H
#define Feta_PTT_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Feta_PTT Declaration
\*---------------------------------------------------------------------------*/

class Feta_PTT
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

            //- Extensibility parameter
            dimensionedScalar epsilon_;

            //- Relaxation time
            dimensionedScalar lambda_;

            //- Chains slip in relation to the continue environment
            dimensionedScalar zeta_;

            //- Parameters of the modified Ellis model
            dimensionedScalar A_;

            dimensionedScalar a_;

            dimensionedScalar b_;


        //- Effective zero shear rate polymer viscosity
        volScalarField etaPEff_;

        //- Relaxation time
        volScalarField lambdaEff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Feta_PTT(const Feta_PTT&);

        //- Disallow default bitwise assignment
        void operator=(const Feta_PTT&);


public:

    //- Runtime type information
    TypeName("Feta-PTT");

    // Constructors

        //- Construct from components
        Feta_PTT
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~Feta_PTT()
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
