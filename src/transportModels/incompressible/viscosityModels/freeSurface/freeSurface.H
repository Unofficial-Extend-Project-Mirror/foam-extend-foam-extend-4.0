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

        virtual ~freeSurface()
        {}


    // Member Functions

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = gamma_*nu1_ + (scalar(1) - gamma_)*nu2_;
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
