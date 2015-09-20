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
    S-MDCPP

Description
    S-MDCPP non linear viscoelastic fluid model:
    Wei Wang, Xikui Li and Xianhong Han,
    Journal of Non-Newtonian Fluid Mechanics, 165,1480-1493, 2010.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    S_MDCPP.C

\*---------------------------------------------------------------------------*/

#ifndef S_MDCP_H
#define S_MDCP_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class S_MDCPP Declaration
\*---------------------------------------------------------------------------*/

class S_MDCPP
:
    public viscoelasticLaw
{
    // Private data

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
            dimensionedScalar zeta_;

            //- Relaxation time of the backbone tube orientation
            dimensionedScalar lambdaOb_;

            //- Relaxation time for the stretch
            dimensionedScalar lambdaOs_;

            //- Amount of arms at the end of a backbone
            dimensionedScalar q_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        S_MDCPP(const S_MDCPP&);

        //- Disallow default bitwise assignment
        void operator=(const S_MDCPP&);


public:

    //- Runtime type information
    TypeName("S_MDCPP");

    // Constructors

        //- Construct from components
        S_MDCPP
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~S_MDCPP()
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
