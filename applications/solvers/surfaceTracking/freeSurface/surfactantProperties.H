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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    surfactantProperties

Description

SourceFiles
    surfactantProperties.H

\*---------------------------------------------------------------------------*/


#ifndef SurfactantProperties_H
#define SurfactantProperties_H


#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class surfactantProperties
{
    // Private data
    
        //- Surfactant concetration in the bulk of fluid
        dimensionedScalar surfactBulkConc_;

        //- Saturated surfactant concentration on the free-surface
        dimensionedScalar surfactSaturatedConc_;

        //- Adsorption coefficient of surfactant
        dimensionedScalar surfactAdsorptionCoeff_;

        //- Desorption coefficient of surfactant
        dimensionedScalar surfactDesorptionCoeff_;

        //- Diffusion coefficient of surfactant in the bulk of fluid
        dimensionedScalar surfactBulkDiffusion_;

        //- Diffusion coefficient of surfactant at the free-surface
        dimensionedScalar surfactDiffusion_;

        //- Temperature of surfactant at the free-surface
        dimensionedScalar surfactT_;

        //- Universal gas constant
        dimensionedScalar surfactR_;

        //- Equilibrium surfactant concentration at the free-surface
        dimensionedScalar surfactEquilibriumConc_;

        //- Is the surfactant soluble
        Switch soluble_;

public:

    // Constructors
    surfactantProperties(const dictionary& dict)
    :
        surfactBulkConc_(dict.lookup("bulkConc")),
        surfactSaturatedConc_(dict.lookup("saturatedConc")),
        surfactAdsorptionCoeff_(dict.lookup("adsorptionCoeff")),
        surfactDesorptionCoeff_(dict.lookup("desorptionCoeff")),
        surfactBulkDiffusion_(dict.lookup("bulkDiffusion")),
        surfactDiffusion_(dict.lookup("diffusion")),
        surfactT_(dict.lookup("temperature")),
        surfactR_("R", dimGasConstant*dimMass/dimMoles, 8.3144),
        surfactEquilibriumConc_
        (
            surfactSaturatedConc_
           /(
                1.0 
              + surfactDesorptionCoeff_
               /surfactBulkConc_
            )
        ),
        soluble_(dict.lookup("soluble"))
    {}


    // Member function

        //- Return surfactant concentration in the bulk of fluid
        const dimensionedScalar& surfactBulkConc() const
        {
            return surfactBulkConc_;
        }

        //- Return saturated surfactant concentration at the free-surface
        const dimensionedScalar& surfactSaturatedConc() const
        {
            return surfactSaturatedConc_;
        }

        //- Return surfactant adsorption coefficient
        const dimensionedScalar& surfactAdsorptionCoeff() const
        {
            return surfactAdsorptionCoeff_;
        }

        //- Return surfactant desorption coefficient
        const dimensionedScalar& surfactDesorptionCoeff() const
        {
            return surfactDesorptionCoeff_;
        }

        //- Return diffusion coefficient of the surfactant in the bulk of fluid
        const dimensionedScalar& surfactBulkDiffusion() const
        {
            return surfactBulkDiffusion_;
        }

        //- Return diffusion coefficient of the surfactant at the free-surface
        const dimensionedScalar& surfactDiffusion() const
        {
            return surfactDiffusion_;
        }

        //- Return surfactant temeprature
        const dimensionedScalar& surfactT() const
        {
            return surfactT_;
        }

        //- Return universal gas constant
        const dimensionedScalar& surfactR() const
        {
            return surfactR_;
        }

        //- Return equilibrium surfactant concentration at the free-surface
        const dimensionedScalar& surfactEquilibriumConc() const
        {
            return surfactEquilibriumConc_;
        }

        //- Is the surfactant soluble
        Switch soluble() const
        {
            return soluble_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
