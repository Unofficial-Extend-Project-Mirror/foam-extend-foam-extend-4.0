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
    plasticityModel

Description
    Class for plasticity model.

SourceFiles
    plasticityModel.C

\*---------------------------------------------------------------------------*/

#ifndef plasticityModel_H
#define plasticityModel_H

#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plasticityModel Declaration
\*---------------------------------------------------------------------------*/

class plasticityModel
:
    public rheologyModel
{
    // Private data

        const volTensorField& gradDU_;

        const volSymmTensorField& epsilon_;

        dictionary plasticityModelCoeffs_;

        volScalarField beta_;

        volScalarField sigmaY_;

        volSymmTensorField DEpsilonP_;

        volScalarField mu_;

        volScalarField lambda_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        plasticityModel(const plasticityModel&);

        //- Disallow default bitwise assignment
        void operator=(const plasticityModel&);


public:

    //- Runtime type information
    TypeName("plasticityModel");


    // Constructors

        //- Construct from components
        plasticityModel
        (
            const volTensorField& gradDU,
            const volSymmTensorField& epsilon,
            const volSymmTensorField& sigma
        );


    // Destructor

        virtual ~plasticityModel();


    // Member Functions

        // Access

        //- Return plastic strain increment
        const volSymmTensorField& DEpsilonP() const
        {
            return DEpsilonP_;
        }
    
        const volScalarField& newMu() const
        {
            return mu_;
        }
    
        const volScalarField& newLambda() const
        {
            return lambda_;
        }
    
        //- Correct plastic strain increment
        virtual void correct();

        //- Update yield stress
        void updateYieldStress();

        //- Read plasticityProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
