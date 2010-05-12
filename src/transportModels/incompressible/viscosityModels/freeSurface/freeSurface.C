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
