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

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(freeSurface, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        freeSurface,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::freeSurface::freeSurface
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    freeSurfaceCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nu1_(freeSurfaceCoeffs_.lookup("nu1")),
    nu2_(freeSurfaceCoeffs_.lookup("nu2")),

    gamma_
    (
        IOobject
        (
            "gamma",
            U_.time().timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh()
    ),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gamma_*nu1_ + (scalar(1) - gamma_)*nu2_
    )
{};


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::freeSurface::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    freeSurfaceCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    freeSurfaceCoeffs_.lookup("nu1") >> nu1_;
    freeSurfaceCoeffs_.lookup("nu2") >> nu2_;

    return true;
}


// ************************************************************************* //
