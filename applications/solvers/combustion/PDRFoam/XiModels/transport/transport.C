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

#include "transport.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(transport, 0);
    addToRunTimeSelectionTable(XiModel, transport, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::transport::transport
(
    const dictionary& XiProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su,
    const volScalarField& rho,
    const volScalarField& b,
    const surfaceScalarField& phi
)
:
    XiModel(XiProperties, thermo, turbulence, Su, rho, b, phi),
    XiShapeCoef(readScalar(XiModelCoeffs_.lookup("XiShapeCoef"))),
    XiEqModel_(XiEqModel::New(XiProperties, thermo, turbulence, Su)),
    XiGModel_(XiGModel::New(XiProperties, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiModels::transport::~transport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::transport::Db() const
{
    return XiGModel_->Db();
}


void Foam::XiModels::transport::correct
(
    const fv::convectionScheme<scalar>& mvConvection
)
{
    volScalarField XiEqEta = XiEqModel_->XiEq();
    volScalarField GEta = XiGModel_->G();

    volScalarField R = GEta*XiEqEta/(XiEqEta - 0.999);

    volScalarField XiEqStar = R/(R - GEta);

    volScalarField XiEq =
        1.0 + (1.0 + (2*XiShapeCoef)*(0.5 - b_))*(XiEqStar - 1.0);

    volScalarField G = R*(XiEq - 1.0)/XiEq;

    const objectRegistry& db = b_.db();
    const volScalarField& betav = db.lookupObject<volScalarField>("betav");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");
    const surfaceScalarField& phiSt =
        db.lookupObject<surfaceScalarField>("phiSt");
    const volScalarField& Db = db.lookupObject<volScalarField>("Db");
    const surfaceScalarField& nf = db.lookupObject<surfaceScalarField>("nf");

    surfaceScalarField phiXi
    (
        "phiXi",
        phiSt
      + (
          - fvc::interpolate(fvc::laplacian(Db, b_)/mgb)*nf
          + fvc::interpolate(rho_)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf
        )
    );

    solve
    (
        betav*fvm::ddt(rho_, Xi_)
      + mvConvection.fvmDiv(phi_, Xi_)
      + fvm::div(phiXi, Xi_)
      - fvm::Sp(fvc::div(phiXi), Xi_)
     ==
        betav*rho_*R
      - fvm::Sp(betav*rho_*(R - G), Xi_)
    );

    // Correct boundedness of Xi
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    Xi_.max(1.0);
    Xi_ = min(Xi_, 2.0*XiEq);
}


bool Foam::XiModels::transport::read(const dictionary& XiProperties)
{
    XiModel::read(XiProperties);

    XiModelCoeffs_.lookup("XiShapeCoef") >> XiShapeCoef;

    return true;
}


// ************************************************************************* //
