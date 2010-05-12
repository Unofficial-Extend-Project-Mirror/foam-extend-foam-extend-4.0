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

#include "LESModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LESModel, 0);
defineRunTimeSelectionTable(LESModel, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void LESModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

LESModel::LESModel
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermoPhysicalModel
)
:
    IOdictionary
    (
        IOobject
        (
            "LESProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),
    phi_(phi),
    thermoPhysicalModel_(thermoPhysicalModel),

    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subDict(type + "Coeffs")),

    k0_("k0", dimVelocity*dimVelocity, SMALL),

    delta_(LESdelta::New("delta", U.mesh(), *this))
{
    if (found("k0"))
    {
        lookup("k0") >> k0_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LESModel::correct(const tmp<volTensorField>&)
{
    delta_().correct();
}


void LESModel::correct()
{
    correct(fvc::grad(U_));
}


bool LESModel::read()
{
    if (regIOobject::read())
    {
        coeffDict_ = subDict(type() + "Coeffs");

        delta_().read(*this);

        if (found("k0"))
        {
            lookup("k0") >> k0_;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
