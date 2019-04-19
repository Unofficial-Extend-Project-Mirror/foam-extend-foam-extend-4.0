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

#include "LienCubicKELowRe.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienCubicKELowRe, 0);
addToRunTimeSelectionTable(RASModel, LienCubicKELowRe, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LienCubicKELowRe::LienCubicKELowRe
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
    CmuWall_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    kappa_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Am_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Am",
            coeffDict_,
            0.016
        )
    ),
    Aepsilon_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Aepsilon",
            coeffDict_,
            0.263
        )
    ),
    Amu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Amu",
            coeffDict_,
            0.00222
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

    y_(mesh_),

    eta_(k_/epsilon_*sqrt(2.0*magSqr(symm(fvc::grad(U_))))),
    ksi_(k_/epsilon_*sqrt(2.0*magSqr(skew(fvc::grad(U_))))),
    Cmu_(2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_))),
    fEta_(A2_ + pow(eta_, 3.0)),

    C5viscosity_
    (
       -2*pow3(Cmu_)*pow4(k_)/pow3(epsilon_)*
       (
           magSqr(twoSymm(fvc::grad(U_)))
         - magSqr(2*skew(fvc::grad(U_)))
        )
    ),

    yStar_(sqrt(k_)*y_/nu() + SMALL),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateLowReNut("nut", mesh_)
    ),

    nonlinearStress_
    (
        "nonlinearStress",
        symm
        (
            // quadratic terms
            pow3(k_)/sqr(epsilon_)
           *(
                Ctau1_/fEta_
               *(
                    (fvc::grad(U_) & fvc::grad(U_))
                  + (fvc::grad(U_) & fvc::grad(U_))().T()
                )
              + Ctau2_/fEta_*(fvc::grad(U_) & T(fvc::grad(U_)))
              + Ctau3_/fEta_*(T(fvc::grad(U_)) & fvc::grad(U_))
            )
            // cubic term C4
          - 20.0*pow(k_, 4.0)/pow3(epsilon_)
           *pow3(Cmu_)
           *(
                ((fvc::grad(U_) & fvc::grad(U_)) & T(fvc::grad(U_)))
              + ((fvc::grad(U_) & T(fvc::grad(U_))) & T(fvc::grad(U_)))
              - ((T(fvc::grad(U_)) & fvc::grad(U_)) & fvc::grad(U_))
              - ((T(fvc::grad(U_)) & T(fvc::grad(U_))) & fvc::grad(U_))
            )
            // cubic term C5, explicit part
          + min
            (
                C5viscosity_,
                dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
            )*fvc::grad(U_)
        )
    )
{
    nut_ = Cmu_
       *(
            scalar(1) - exp(-Am_*yStar_))
           /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL
        )
       *sqr(k_)/(epsilon_ + epsilonSmall_)
        // cubic term C5, implicit part
      + max
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        );

    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LienCubicKELowRe::R() const
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


tmp<volSymmTensorField> LienCubicKELowRe::devReff() const
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


tmp<fvVectorMatrix> LienCubicKELowRe::divDevReff() const
{
    const volScalarField nuEffective = nuEff();

    return
    (
        fvc::div(nonlinearStress_)
      - fvm::laplacian(nuEffective, U_)
      - (fvc::grad(U_) & fvc::grad(nuEffective))
    );
}


bool LienCubicKELowRe::read()
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
        CmuWall_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Am_.readIfPresent(coeffDict());
        Aepsilon_.readIfPresent(coeffDict());
        Amu_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void LienCubicKELowRe::correct()
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

    if (mesh_.changing())
    {
        y_.correct();
    }

    // Changed return type for gradient cacheing.  HJ, 22/Apr/2016
    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    // generation term
    volScalarField S2 = symm(gradU) && gradU;

    yStar_ = sqrt(k_)*y_/nu() + SMALL;
    volScalarField Rt = sqr(k_)/(nu()*epsilon_);

    volScalarField fMu =
        (scalar(1) - exp(-Am_*yStar_))
       /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL);

    volScalarField f2 = scalar(1) - 0.3*exp(-sqr(Rt));

    volScalarField G =
        Cmu_*fMu*sqr(k_)/epsilon_*S2 - (nonlinearStress_ && gradU);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
        // E-term
      + C2_*f2*pow(Cmu_, 0.75)*pow(k_, scalar(0.5))
       /(kappa_*y_*(scalar(1) - exp(-Aepsilon_*yStar_)))
       *exp(-Amu_*sqr(yStar_))*epsilon_
      - fvm::Sp(C2_*f2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "LienCubicKELowReSetWallDissipation.H"
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
    fEta_ = A2_ + pow3(eta_);

    C5viscosity_ =
      - 2.0*pow3(Cmu_)*pow(k_, 4.0)/pow3(epsilon_)
       *(magSqr(gradU + gradU.T()) - magSqr(gradU - gradU.T()));

    nut_ =
        Cmu_*fMu*sqr(k_)/epsilon_
        // C5 term, implicit
      + max
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        );
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    nonlinearStress_ = symm
    (
        // quadratic terms
        pow(k_, 3.0)/sqr(epsilon_)
       *(
            Ctau1_/fEta_
           *(
                (gradU & gradU)
              + (gradU & gradU)().T()
            )
          + Ctau2_/fEta_*(gradU & gradU.T())
          + Ctau3_/fEta_*(gradU.T() & gradU)
        )
        // cubic term C4
      - 20.0*pow(k_, 4.0)/pow(epsilon_, 3.0)
       *pow(Cmu_, 3.0)
       *(
            ((gradU & gradU) & gradU.T())
          + ((gradU & gradU.T()) & gradU.T())
          - ((gradU.T() & gradU) & gradU)
          - ((gradU.T() & gradU.T()) & gradU)
        )
        // cubic term C5, explicit part
      + min
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        )*gradU
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
