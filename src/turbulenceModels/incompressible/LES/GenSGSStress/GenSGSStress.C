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

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GenSGSStress::GenSGSStress
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(word("GenSGSStress"), U, phi, transport),

    ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            coeffDict_,
            1.048
        )
    ),

    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            coeffDict_,
            0.0
        )
    ),

    B_
    (
        IOobject
        (
            "B",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nu(),
        B_.boundaryField().types()
    )
{
    if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
    {
        FatalErrorIn
        (
            "GenSGSStress::GenSGSStress"
            "(const volVectorField& U, const surfaceScalarField& phi,"
            "transportModel& lamTransportModel)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> GenSGSStress::devBeff() const
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
            B_ - nu()*dev(twoSymm(fvc::grad(U())))
        )
    );
}


tmp<fvVectorMatrix> GenSGSStress::divDevBeff
(
    volVectorField& U
) const
{
    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(B_ + couplingFactor_*nuSgs_*fvc::grad(U))
          + fvc::laplacian
            (
                (1.0 - couplingFactor_)*nuSgs_, U, "laplacian(nuEff,U)"
            )
          - fvm::laplacian(nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::div(B_)
          + fvc::laplacian(nuSgs_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
}


bool GenSGSStress::read()
{
    if (LESModel::read())
    {
        ce_.readIfPresent(coeffDict());

        couplingFactor_.readIfPresent(coeffDict());

        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorIn("GenSGSStress::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1"
                << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
