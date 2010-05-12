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

#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radiationModel, 0);
        defineRunTimeSelectionTable(radiationModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.mesh().time().constant(),
            T.mesh().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    T_(T),
    mesh_(T.mesh()),
    radiation_(lookup("radiation")),
    radiationModelCoeffs_(subDict(type + "Coeffs")),
    absorptionEmission_(absorptionEmissionModel::New(*this, mesh_)),
    scatter_(scatterModel::New(*this, mesh_))
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::radiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("radiation") >> radiation_;
        radiationModelCoeffs_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::Sh
(
    basicThermo& thermo
) const
{
    volScalarField& h = thermo.h();
    const volScalarField cp = thermo.Cp();
    const volScalarField T3 = pow3(T_);

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/cp, h)
      - Rp()*T3*(T_ - 4.0*h/cp)
    );
}


// ************************************************************************* //
