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

#include "orthotropicLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformField.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(orthotropicLinearElastic, 0);
    addToRunTimeSelectionTable
    (rheologyLaw, orthotropicLinearElastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::orthotropicLinearElastic::orthotropicLinearElastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    rho_(dict.lookup("rho")),
    Ex_(dict.lookup("Ex")),
    Ey_(dict.lookup("Ey")),
    Ez_(dict.lookup("Ez")),
    E_(Ex_),
    nuxy_(dict.lookup("nuxy")),
    nuyz_(dict.lookup("nuyz")),
    nuzx_(dict.lookup("nuzx")),
    nuyx_(nuxy_*Ey_/Ex_),
    nuzy_(nuyz_*Ez_/Ey_),
    nuxz_(nuzx_*Ex_/Ez_),
    nu_(nuxy_),
    Gxy_(dict.lookup("Gxy")),
    Gyz_(dict.lookup("Gyz")),
    Gzx_(dict.lookup("Gzx")),
    J_(
        (1 - nuxy_*nuyx_ - nuyz_*nuzy_ - nuzx_*nuxz_ - 2*nuyx_*nuzy_*nuxz_)
        / (Ex_*Ey_*Ez_)
        ),
    A11_( (1 - nuyz_*nuzy_)/(J_*Ey_*Ez_) ),
    A22_( (1 - nuxz_*nuzx_)/(J_*Ex_*Ez_) ),
    A33_( (1 - nuyx_*nuxy_)/(J_*Ey_*Ex_) ),
    A12_( (nuxy_ + nuzy_*nuxz_)/(J_*Ex_*Ez_) ),
    A31_( (nuzx_ + nuyx_*nuzy_)/(J_*Ey_*Ez_) ),
    A23_( (nuyz_ + nuyx_*nuxz_)/(J_*Ex_*Ey_) ),
    A44_( 2*Gxy_ ),
    A55_( 2*Gyz_ ),
    A66_( 2*Gzx_ ),
    C_(
       IOobject
       (
           "rheologyLawStoredC",
           mesh().time().timeName(),
           mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE
           ),
       mesh(),
       dimensionedSymmTensor4thOrder
       (
           "zero",
           dimForce/dimArea,
           symmTensor4thOrder(A11_.value(), A12_.value(), A31_.value(),
                              A22_.value(), A23_.value(),
                              A33_.value(),
                              A44_.value(),
                              A55_.value(),
                              A66_.value())
           ),
       zeroGradientFvPatchSymmTensor4thOrderField::typeName
        ),
    matDir_(
        IOobject
        (
            "materialDirections",
         mesh().time().timeName(),
         mesh(),
         IOobject::MUST_READ,
         IOobject::NO_WRITE
         ),
        mesh()
        )
{
  //- check material properties lie within physical constraints
  //- ref Abaqus analysis user's manual orthotropic material
    Info<< "\tChecking physical constraints on the orthotropic"
        << " material properties" << endl;

  //- E and G should be greater than zero
  if (Ex_.value() < 0.0 || Ey_.value() < 0.0 || Ez_.value() < 0.0
     || Gxy_.value() < 0.0 || Gyz_.value() < 0.0 || Gzx_.value() < 0.0)
    {
      FatalError << "Ex, Ey, Ez, Gxy, Gyz, Gzx should all be greater than zero!"
         << exit(FatalError);
    }

  //- restriction on Poisson's ratio
  if (mag(nuxy_.value()) >= sqrt(Ex_.value()/Ey_.value())
     || mag(nuyz_.value()) >= sqrt(Ey_.value()/Ez_.value())
     || mag(nuzx_.value()) >= sqrt(Ez_.value()/Ex_.value()))
    {
      FatalError << "mag(nuij) should be less sqrt(Ei/Ej)"
         << exit(FatalError);
    }

  if (
      dimensionedScalar
      (1 - nuxy_*nuyx_ - nuyz_*nuzy_
       - nuzx_*nuxz_ - 2*nuyx_*nuzy_*nuxz_).value()
      <= 0 )
    {
      FatalError
          << "(1 - nuxy*nuyx - nuyz*nuzy "
          << "- nuzx*nuxz - 2*nuyx*nuzy*nuxz) should be greater than zero!"
         << exit(FatalError);
    }


  Info<< "\tRotating local material properties to global coordinate system"
      << endl;
  //- rotate tensors
  volTensorField R
    (
     IOobject
     (
      "R",
      mesh().time().timeName(),
      mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh(),
     dimensionedTensor("zero", dimless, tensor::zero),
     zeroGradientFvPatchTensorField::typeName
     );

  //- make sure matDir_ are unit directions
  forAll(matDir_, celli)
    {
      {
    scalar magVec =
      mag(vector(matDir_[celli][0], matDir_[celli][1], matDir_[celli][2]));
    matDir_[celli][0] /= magVec;
    matDir_[celli][1] /= magVec;
    matDir_[celli][2] /= magVec;
      }
      {
    scalar magVec =
      mag(vector(matDir_[celli][3], matDir_[celli][4], matDir_[celli][5]));
    matDir_[celli][3] /= magVec;
    matDir_[celli][4] /= magVec;
    matDir_[celli][5] /= magVec;
      }
      {
    scalar magVec =
      mag(vector(matDir_[celli][6], matDir_[celli][7], matDir_[celli][8]));
    matDir_[celli][6] /= magVec;
    matDir_[celli][7] /= magVec;
    matDir_[celli][8] /= magVec;
      }
    }

  //- global axes
  vector e0(1,0,0);
  vector e1(0,1,0);
  vector e2(0,0,1);

  forAll(R, celli)
    {
      // R_ij = xold_i & xnew_i;
      {
    vector mD(matDir_[celli][0], matDir_[celli][1], matDir_[celli][2]);
    R[celli][0] = e0 & mD;
    R[celli][1] = e1 & mD;
    R[celli][2] = e2 & mD;
      }
      {
    vector mD(matDir_[celli][3], matDir_[celli][4], matDir_[celli][5]);
    R[celli][3] = e0 & mD;
    R[celli][4] = e1 & mD;
    R[celli][5] = e2 & mD;
      }
      {
    vector mD(matDir_[celli][6], matDir_[celli][7], matDir_[celli][8]);
    R[celli][6] = e0 & mD;
    R[celli][7] = e1 & mD;
    R[celli][8] = e2 & mD;
      }
    }

  //Info << "R is " << R.internalField() << endl;
  R.correctBoundaryConditions();
  //R.write();

  //- rotate C to global corrdinate system
  //C_.correctBoundaryConditions();
  C_ = transform(R, C_);
  C_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::orthotropicLinearElastic::~orthotropicLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::orthotropicLinearElastic::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::orthotropicLinearElastic::E() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            E_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::orthotropicLinearElastic::nu() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}

Foam::tmp<Foam::volScalarField> Foam::orthotropicLinearElastic::Ep() const
{
//     notImplemented(type() + "::Ep()");

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Ep",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroEp", dimForce/dimArea, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::orthotropicLinearElastic::sigmaY() const
{
//     notImplemented(type() + "::sigmaY()");

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaY",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroSigmaY", dimForce/dimArea, GREAT),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Foam::scalar Foam::orthotropicLinearElastic::sigmaY
(const scalar epsilonPEq, const label cellID) const
{
  return GREAT;
}

//Foam::tmp<Foam::volTensorField> Foam::orthotropicLinearElastic::K() const
Foam::tmp<Foam::volDiagTensorField> Foam::orthotropicLinearElastic::K() const
{
  tmp<volDiagTensorField> tresult
    (
        new volDiagTensorField
        (
            IOobject
            (
                "K",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedDiagTensor("zero", A11_.dimensions(), diagTensor::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

  volDiagTensorField& K = tresult();

  forAll(K, celli)
    {
      K[celli].xx() = C_[celli].xxxx();
      K[celli].yy() = C_[celli].yyyy();
      K[celli].zz() = C_[celli].zzzz();
    }

  K.correctBoundaryConditions();

  return tresult;
}

Foam::tmp<Foam::volSymmTensor4thOrderField>
Foam::orthotropicLinearElastic::C() const
{
  tmp<volSymmTensor4thOrderField> tresult
    (
        new volSymmTensor4thOrderField
        (
     C_
     )
    );

  volSymmTensor4thOrderField& result = tresult();

  result.correctBoundaryConditions();

  return tresult;
}
// ************************************************************************* //
