/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2008 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    rheologyModel

Description
    Material rheology for solids.

\*---------------------------------------------------------------------------*/

#include "rheologyModel.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rheologyModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rheologyModel::rheologyModel
(
    const volSymmTensorField& sigma
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
    planeStress_(lookup("planeStress")),
    lawPtr_(rheologyLaw::New("law", sigma_, subDict("rheology")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Return first Lame's coefficient
tmp<volScalarField> rheologyModel::mu() const
{
    volScalarField lawE = lawPtr_->E();
    volScalarField lawNu = lawPtr_->nu();

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
tmp<volScalarField> rheologyModel::lambda() const
{
    volScalarField lawE = lawPtr_->E();
    volScalarField lawNu = lawPtr_->nu();

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
tmp<volScalarField> rheologyModel::threeK() const
{
    volScalarField lawRho = lawPtr_->rho();
    volScalarField lawE = lawPtr_->E();
    volScalarField lawNu = lawPtr_->nu();

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
tmp<volScalarField> rheologyModel::mu(scalar t) const
{
    volScalarField lawE = lawPtr_->E(t);
    volScalarField lawNu = lawPtr_->nu(t);

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
tmp<volScalarField> rheologyModel::lambda(scalar t) const
{
    volScalarField lawE = lawPtr_->E(t);
    volScalarField lawNu = lawPtr_->nu(t);

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


void rheologyModel::correct()
{
    lawPtr_->correct();
}


bool rheologyModel::read()
{
    if (regIOobject::read())
    {
        lookup("planeStress") >> planeStress_;
        lawPtr_ = rheologyLaw::New("law", sigma_, subDict("rheology"));

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
