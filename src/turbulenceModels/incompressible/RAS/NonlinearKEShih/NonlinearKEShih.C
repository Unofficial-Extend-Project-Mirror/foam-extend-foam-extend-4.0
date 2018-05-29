/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "NonlinearKEShih.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NonlinearKEShih, 0);
addToRunTimeSelectionTable(RASModel, NonlinearKEShih, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NonlinearKEShih::NonlinearKEShih
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    C1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmak_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    A1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "A1",
            coeffDict_,
            1.25
        )
    ),
    A2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "A2",
            coeffDict_,
            1000.0
        )
    ),
    Ctau1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ctau1",
            coeffDict_,
            -4.0
        )
    ),
    Ctau2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ctau2",
            coeffDict_,
            13.0
        )
    ),
    Ctau3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ctau3",
            coeffDict_,
            -2.0
        )
    ),
    alphaKsi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaKsi",
            coeffDict_,
            0.9
        )
    ),

    kappa_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa_",
            coeffDict_,
            0.41
        )
    ),
    E_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "E",
            coeffDict_,
            9.8
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    eta_(k_/epsilon_*sqrt(2.0*magSqr(symm(fvc::grad(U_))))),
    ksi_(k_/epsilon_*sqrt(2.0*magSqr(skew(fvc::grad(U_))))),
    Cmu_(2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_))),
    fEta_(A2_ + pow(eta_, 3.0)),

    nut_("nut", Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_)),

    nonlinearStress_
    (
        "nonlinearStress",
        symm
        (
            pow(k_, 3.0)/sqr(epsilon_)
           *(
                Ctau1_/fEta_
               *(
                    (fvc::grad(U_) & fvc::grad(U_))
                  + T(fvc::grad(U_) & fvc::grad(U_))
                )
              + Ctau2_/fEta_*(fvc::grad(U_) & T(fvc::grad(U_)))
              + Ctau3_/fEta_*(T(fvc::grad(U_)) & fvc::grad(U_))
            )
        )
    )
{
#   include "wallNonlinearViscosityI.H"

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> NonlinearKEShih::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + nonlinearStress_,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> NonlinearKEShih::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_))) + nonlinearStress_
        )
    );
}


tmp<fvVectorMatrix> NonlinearKEShih::divDevReff() const
{
    const volScalarField nuEffective = nuEff();

    return
    (
        fvc::div(nonlinearStress_)
      - fvm::laplacian(nuEffective, U_)
      - (fvc::grad(U_) & fvc::grad(nuEffective))
    );
}


bool NonlinearKEShih::read()
{
    if (RASModel::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        A1_.readIfPresent(coeffDict());
        A2_.readIfPresent(coeffDict());
        Ctau1_.readIfPresent(coeffDict());
        Ctau2_.readIfPresent(coeffDict());
        Ctau3_.readIfPresent(coeffDict());
        alphaKsi_.readIfPresent(coeffDict());

        kappa_.readIfPresent(coeffDict());
        E_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void NonlinearKEShih::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(epsilon_, epsilon0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    // Changed return type for gradient cacheing.  HJ, 22/Apr/2016
    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    // generation term
    volScalarField S2 = symm(gradU) && gradU;

    volScalarField G
    (
        "RASModel::G",
        Cmu_*sqr(k_)/epsilon_*S2
      - (nonlinearStress_ && gradU)
    );

#   include "nonLinearWallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity

    eta_ = k_/epsilon_*sqrt(2.0*magSqr(symm(gradU)));
    ksi_ = k_/epsilon_*sqrt(2.0*magSqr(skew(gradU)));
    Cmu_ = 2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_));
    fEta_ = A2_ + pow(eta_, 3.0);

    nut_ = Cmu_*sqr(k_)/epsilon_;

#   include "wallNonlinearViscosityI.H"

    nonlinearStress_ = symm
    (
        pow(k_, 3.0)/sqr(epsilon_)
       *(
            Ctau1_/fEta_
           *(
                (gradU & gradU)
              + T(gradU & gradU)
            )
          + Ctau2_/fEta_*(gradU & gradU.T())
          + Ctau3_/fEta_*(gradU.T() & gradU)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
