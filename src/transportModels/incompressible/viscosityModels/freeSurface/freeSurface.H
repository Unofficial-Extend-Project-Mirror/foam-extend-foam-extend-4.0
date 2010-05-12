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
    Foam::viscosityModels::freeSurface

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    freeSurface.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurface_H
#define freeSurface_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class freeSurface
:
    public viscosityModel
{
    // Private data

        dictionary freeSurfaceCoeffs_;

        //- Viscosity of two materials
        dimensionedScalar nu1_;
        dimensionedScalar nu2_;

        //- Free surface indicator variable
        volScalarField gamma_;

        //- Viscosity
        volScalarField nu_;

public:

    //- Runtime type information
    TypeName("freeSurface");


    // Constructors

        //- Construct from components
        freeSurface
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        ~freeSurface()
        {}


    // Member Functions

        //- Return the laminar viscosity
        const volScalarField& nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = gamma_*nu1_ + (scalar(1) - gamma_)*nu2_;
	}

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
