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

#include "SpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmaras, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::fv1() const
{
    volScalarField chi3 = pow3(nuTilda_/nu());
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SpalartAllmaras::fv2() const
{
    volScalarField chi = nuTilda_/nu();
    //return scalar(1) - chi/(scalar(1) + chi*fv1());
    return 1.0/pow3(scalar(1) + chi/Cv2_);
}


tmp<volScalarField> SpalartAllmaras::fv3() const
{
    volScalarField chi = nuTilda_/nu();
    volScalarField chiByCv2 = (1/Cv2_)*chi;

    return
        (scalar(1) + chi*fv1())
       *(1/Cv2_)
       *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
       /pow3(scalar(1) + chiByCv2);
}


tmp<volScalarField> SpalartAllmaras::fw(const volScalarField& Stilda) const
{
    volScalarField r = min
    (
        nuTilda_
       /(
           max(Stilda, dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
          *sqr(kappa_*dTilda_)
        ),
        scalar(10.0)
    );
    r.boundaryField() == 0.0;

    volScalarField g = r + Cw2_*(pow6(r) - r);

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(typeName, U, phi, transport),


    alphaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaNut",
            coeffDict(),
            1.5
        )
    ),
    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            coeffDict(),
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            coeffDict(),
            0.622
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            coeffDict(),
            7.1
        )
    ),
    Cv2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv2",
            coeffDict(),
            5.0
        )
    ),
    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            coeffDict(),
            0.65
        )
    ),
    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict(),
            0.07
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            *this,
            0.4187
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + alphaNut_*(1.0 + Cb2_)),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict(),
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict(),
            2.0
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    dTilda_(min(CDES_*delta(), wallDist(mesh_).y())),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SpalartAllmaras::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    if (mesh_.changing())
    {
        dTilda_ = min(CDES_*delta(), wallDist(mesh_).y());
    }

    volScalarField Stilda =
        fv3()*::sqrt(2.0)*mag(skew(gradU)) + fv2()*nuTilda_/sqr(kappa_*dTilda_);

    solve
    (
        fvm::ddt(nuTilda_)
      + fvm::div(phi(), nuTilda_)
      - fvm::laplacian
        (
            alphaNut_*(nuTilda_ + nu()),
            nuTilda_,
            "laplacian(DnuTildaEff,nuTilda)"
        )
      - alphaNut_*Cb2_*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*Stilda*nuTilda_
      - fvm::Sp(Cw1_*fw(Stilda)*nuTilda_/sqr(dTilda_), nuTilda_)
    );

    bound(nuTilda_, dimensionedScalar("zero", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    nuSgs_.internalField() = fv1()*nuTilda_.internalField();
    nuSgs_.correctBoundaryConditions();
}


tmp<volScalarField> SpalartAllmaras::epsilon() const
{
    return 2*nuEff()*magSqr(symm(fvc::grad(U())));
}


tmp<volSymmTensorField> SpalartAllmaras::B() const
{
    return ((2.0/3.0)*I)*k() - nuSgs()*twoSymm(fvc::grad(U()));
}


tmp<volSymmTensorField> SpalartAllmaras::devBeff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> SpalartAllmaras::divDevBeff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U) - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool SpalartAllmaras::read()
{
    if (LESModel::read())
    {
        alphaNut_.readIfPresent(coeffDict());
        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + alphaNut_*(1.0 + Cb2_);
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cv2_.readIfPresent(coeffDict());
        CDES_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());
        kappa_.readIfPresent(*this);

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
