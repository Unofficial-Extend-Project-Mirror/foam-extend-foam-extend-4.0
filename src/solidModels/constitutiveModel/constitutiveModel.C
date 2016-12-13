/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "constitutiveModel.H"
#include "volFields.H"
#include "fvc.H"
#include "solidInterface.H"
#include "solidCohesiveFvPatchVectorField.H"
#include "solidCohesiveFixedModeMixFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(constitutiveModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constitutiveModel::constitutiveModel
(
 const volSymmTensorField& sigma,
 const volVectorField& U
)
:
    IOdictionary
    (
        IOobject
        (
            "rheologyProperties",
            sigma.time().constant(),
            sigma.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sigma_(sigma),
    rheologyLawPtr_(rheologyLaw::New("law", sigma_, subDict("rheology"))),
    cohesiveDictPtr_(NULL),
    cohesiveLawPtr_(NULL),
    planeStress_(lookup("planeStress")),
    solidInterfacePtr_(NULL),
    solidInterfaceActive_(false),
    plasticityStressReturnPtr_(NULL)
{
    Info << "Creating constitutive model" << endl;

    if (rheologyLawPtr_->plasticityModelNeeded())
    {
        plasticityStressReturnPtr_ =
            plasticityStressReturn::New(lookup("plasticityModel"), *this);
    }
    else if (found("plasticityModel"))
    {
        Warning
            << "plasticityModel is specified in rheologyProperties"
            << " but plasticity is not active in the current model"
                << endl;
    }

    // create solidInterface if the case is multiMaterial
    if (rheologyLawPtr_->type() == "multiMaterial")
    {
        // we should figure it out automatically
        // for now we will ask the user
        //word solIntType("smallStrain");
        solidInterfaceActive_ = true;

        word solIntType
        (
            sigma.mesh().solutionDict().subDict
            (
                "solidMechanics"
            ).lookup("solidInterfaceMethod")
        );

        solidInterfacePtr_ =
            new IOReferencer<solidInterface>
            (
                IOobject
                (
                    "solidInterface",
                    sigma.time().timeName(),
                    sigma.db(),
                    IOobject::NO_READ,  /*must be NO_READ*/
                    IOobject::NO_WRITE  /*must be NO_WRITE*/
                ),
                solidInterface::New(solIntType, sigma.mesh(), *this).ptr()
            );
    }

    // create cohesiveLaw if a solidCohesive patch
    forAll(U.boundaryField(), patchi)
    {
        if
        (
            isA<solidCohesiveFvPatchVectorField>(U.boundaryField()[patchi])
         || isA<solidCohesiveFixedModeMixFvPatchVectorField>
            (U.boundaryField()[patchi])
        )
        {
            Info<< "Reading cohesiveProperties because "
                << "a solidCohesive patch exists" << endl;
            cohesiveDictPtr_ =
                new IOdictionary
                (
                    IOobject
                    (
                        "cohesiveProperties",
                        sigma.time().constant(),
                        sigma.db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

            cohesiveLawPtr_ =
                cohesiveLaw::New
                (
                    "law",
                    sigma_,
                    cohesiveDictPtr_->subDict("cohesive")
                );
            break;
        }
    }
}


constitutiveModel::~constitutiveModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void constitutiveModel::correct()
{
    if (plasticityStressReturnPtr_.valid())
    {
        plasticityStressReturnPtr_->correct();
    }
    else
    {
        rheologyLawPtr_->correct();
    }
}

void constitutiveModel::updateYieldStress()
{
    if (plasticityStressReturnPtr_.valid())
    {
        plasticityStressReturnPtr_->updateYieldStress();
    }
}


const volSymmTensorField& constitutiveModel::DEpsilonP() const
{
  return plasticityStressReturnPtr_->DEpsilonP();
}

// Return first Lame's coefficient
tmp<volScalarField> constitutiveModel::mu() const
{
    volScalarField lawE = rheologyLawPtr_->E();
    volScalarField lawNu = rheologyLawPtr_->nu();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            lawE/(2.0*(1.0 + lawNu))
        )
    );
}

// Return second Lame's coefficient
tmp<volScalarField> constitutiveModel::lambda() const
{
    volScalarField lawE = rheologyLawPtr_->E();
    volScalarField lawNu = rheologyLawPtr_->nu();

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - 2.0*lawNu))
            )
        );
    }
}


// Return threeK
tmp<volScalarField> constitutiveModel::threeK() const
{
    volScalarField lawRho = rheologyLawPtr_->rho();
    volScalarField lawE = rheologyLawPtr_->E();
    volScalarField lawNu = rheologyLawPtr_->nu();

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "threeK",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawE/(lawRho*(1 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "threeK",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawE/(lawRho*(1 - 2*lawNu))
            )
        );
    }
}


// Return first Lame's coefficient
tmp<volScalarField> constitutiveModel::mu(scalar t) const
{
    volScalarField lawE = rheologyLawPtr_->E(t);
    volScalarField lawNu = rheologyLawPtr_->nu(t);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            lawE/(2.0*(1.0 + lawNu))
        )
    );
}


// Return first Lame's coefficient
tmp<volScalarField>
constitutiveModel::mu(const volScalarField& epsilonEq) const
{
    volScalarField lawE = rheologyLawPtr_->E(epsilonEq);
    volScalarField lawNu = rheologyLawPtr_->nu();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            lawE/(2.0*(1.0 + lawNu))
        )
    );
}


// Return second Lame's coefficient
tmp<volScalarField> constitutiveModel::lambda(scalar t) const
{
    volScalarField lawE = rheologyLawPtr_->E(t);
    volScalarField lawNu = rheologyLawPtr_->nu(t);

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - 2.0*lawNu))
            )
        );
    }
}


// Return second Lame's coefficient
tmp<volScalarField> constitutiveModel::lambda
(
    const volScalarField& epsilonEq
) const
{
    volScalarField lawE = rheologyLawPtr_->E(epsilonEq);
    volScalarField lawNu = rheologyLawPtr_->nu();

    if (planeStress())
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - lawNu))
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "lambda",
                    sigma_.time().timeName(),
                    sigma_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawNu*lawE/((1.0 + lawNu)*(1.0 - 2.0*lawNu))
            )
        );
    }
}

tmp<volDiagTensorField> constitutiveModel::K() const
{
    return rheologyLawPtr_->K();
}

tmp<volSymmTensor4thOrderField> constitutiveModel::C() const
{
    return rheologyLawPtr_->C();
}

tmp<surfaceScalarField> constitutiveModel::muf() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
                "muf",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(mu(), "mu")
        )
    );
    surfaceScalarField& muf = tresult();

    if (solidInterfaceActive_)
    {
        solInterface().modifyProperties(muf);
    }

    return tresult;
}

tmp<surfaceScalarField> constitutiveModel::lambdaf() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
                "lambdaf",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(lambda(), "lambda")
        )
    );
    surfaceScalarField& lambdaf = tresult();

    if (solidInterfaceActive_)
    {
        solInterface().modifyProperties(lambdaf);
    }

    return tresult;
}

tmp<surfaceScalarField> constitutiveModel::threeKf() const
{
    tmp<surfaceScalarField> tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
                "threeKf",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(threeK(), "threeK")
        )
    );
    surfaceScalarField& threeKf = tresult();

    if (solidInterfaceActive_)
    {
        solInterface().modifyProperties(threeKf);
    }

    return tresult;
}

tmp<surfaceDiagTensorField> constitutiveModel::Kf() const
{
    tmp<surfaceDiagTensorField> tresult
    (
        new surfaceDiagTensorField
        (
            IOobject
            (
                "Kf",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(K(), "K")
        )
    );
    surfaceDiagTensorField& Kf = tresult();

    if (solidInterfaceActive_)
    {
        solInterface().modifyProperties(Kf);
    }

    return tresult;
}

tmp<surfaceSymmTensor4thOrderField> constitutiveModel::Cf() const
{
    tmp<surfaceSymmTensor4thOrderField> tresult
    (
        new surfaceSymmTensor4thOrderField
        (
            IOobject
            (
                "Cf",
                sigma_.time().timeName(),
                sigma_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(C(), "C")
        )
    );
    surfaceSymmTensor4thOrderField& Cf = tresult();

    if (solidInterfaceActive_)
    {
        solInterface().modifyProperties(Cf);
    }

    return tresult;
}


bool constitutiveModel::read()
{
    if (regIOobject::read() && cohesiveDictPtr_->regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
