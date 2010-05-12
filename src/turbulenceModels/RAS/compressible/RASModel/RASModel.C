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

#include "RASModel.H"
#include "wallDist.H"
#include "wallFvPatch.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASModel, 0);
defineRunTimeSelectionTable(RASModel, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void RASModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASModel::RASModel
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    basicThermo& thermophysicalModel
)
:
    IOdictionary
    (
        IOobject
        (
            "RASProperties",
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
    thermophysicalModel_(thermophysicalModel),

    turbulence_(lookup("turbulence")),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subDict(type + "Coeffs")),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            subDict("wallFunctionCoeffs"),
            0.4187
        )
    ),
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "E",
            subDict("wallFunctionCoeffs"),
            9.0
        )
    ),

    yPlusLam_(yPlusLam(kappa_.value(), E_.value())),

    k0_("k0", dimVelocity*dimVelocity, SMALL),
    epsilon0_("epsilon", k0_.dimensions()/dimTime, SMALL),
    epsilonSmall_("epsilonSmall", epsilon0_.dimensions(), SMALL),

    y_(mesh_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar RASModel::yPlusLam(const scalar kappa, const scalar E)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<scalarField> RASModel::yPlus(const label patchNo) const
{
    const fvPatch& curPatch = mesh_.boundary()[patchNo];

    tmp<scalarField> tYp(new scalarField(curPatch.size()));
    scalarField& Yp = tYp();

    if (typeid(curPatch) == typeid(wallFvPatch))
    {
        scalar Cmu(readScalar(coeffDict_.lookup("Cmu")));

        Yp = pow(Cmu, 0.25)
            *y_[patchNo]
            *sqrt(k()().boundaryField()[patchNo].patchInternalField())
           /(
                mu().boundaryField()[patchNo].patchInternalField()
               /rho_.boundaryField()[patchNo]
            );
    }
    else
    {
        WarningIn
        (
            "tmp<scalarField> RASModel::yPlus(const label patchNo) const"
        )   << "Patch " << patchNo << " is not a wall.  Returning blank field"
            << endl;

        Yp.setSize(0);
    }

    return tYp;
}


void RASModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


bool RASModel::read()
{
    if (regIOobject::read())
    {
        lookup("turbulence") >> turbulence_;
        coeffDict_ = subDict(type() + "Coeffs");

        kappa_.readIfPresent(subDict("wallFunctionCoeffs"));
        E_.readIfPresent(subDict("wallFunctionCoeffs"));

        yPlusLam_ = yPlusLam(kappa_.value(), E_.value());

        k0_.readIfPresent(*this);
        epsilon0_.readIfPresent(*this);
        epsilonSmall_.readIfPresent(*this);

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
