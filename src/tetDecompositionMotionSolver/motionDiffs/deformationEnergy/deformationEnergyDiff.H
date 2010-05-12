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
    deformationEnergyDiff

Description
    Deformation energy based motion diffusion.

SourceFiles
    deformationEnergyDiff.C

\*---------------------------------------------------------------------------*/

#ifndef deformationEnergyDiff_H
#define deformationEnergyDiff_H

#include "motionDiff.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class deformationEnergyDiff Declaration
\*---------------------------------------------------------------------------*/

class deformationEnergyDiff
:
    public motionDiff
{
private:

    // Private data
    
        //- Value of exponent
        label exponent_;

        //- NonOrthogonality filed
//         elementScalarField deformationEnergy_;

        //- Diffusion field
        elementScalarField motionGamma_;


   // Private Member Functions

        //- Disallow default bitwise copy construct
        deformationEnergyDiff(const deformationEnergyDiff&);

        //- Disallow default bitwise assignment
        void operator=(const deformationEnergyDiff&);

public:

    //- Runtime type information
    TypeName("deformationEnergy");


    // Constructors

        //- Construct from components
        deformationEnergyDiff
        (
            const tetDecompositionMotionSolver& mSolver
        );


    // Destructor

        virtual ~deformationEnergyDiff();


    // Member Functions

        //- Return diffusivity field
        virtual const elementScalarField& motionGamma() const
        {
            return motionGamma_;
        }

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
