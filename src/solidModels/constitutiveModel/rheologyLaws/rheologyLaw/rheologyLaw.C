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

Class
    rheologyLaw

Description
    Material rheology for solids.

\*---------------------------------------------------------------------------*/

#include "rheologyLaw.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rheologyLaw, 0);
defineRunTimeSelectionTable(rheologyLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rheologyLaw::rheologyLaw
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    name_(name),
    sigma_(sigma)
{}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::volDiagTensorField> Foam::rheologyLaw::K() const
{
  volScalarField J =
      ( (1 - nu()*nu() - nu()*nu() - nu()*nu() - 2*nu()*nu()*nu())
        /(E()*E()*E()) );
  volScalarField A11 = ( (1 - nu()*nu())/(J*E()*E()) );
  volScalarField A22 = ( (1 - nu()*nu())/(J*E()*E()) );
  volScalarField A33 = ( (1 - nu()*nu())/(J*E()*E()) );

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
        dimensionedDiagTensor("K", A11.dimensions(), diagTensor::zero)
        )
    );

  volDiagTensorField& result = tresult();

  result.replace(diagTensor::XX, A11);
  result.replace(diagTensor::YY, A22);
  result.replace(diagTensor::ZZ, A33);

  tresult().correctBoundaryConditions();

  return tresult;
}


Foam::tmp<Foam::volSymmTensor4thOrderField> Foam::rheologyLaw::C() const
{
  volScalarField twoMu = 2.0*E()/(2.0*(1+nu()));
  volScalarField lambda = nu()*E()/((1.0 + nu())*(1.0 - 2.0*nu()));
  volScalarField twoMuLambda = twoMu + lambda;

  // symmTensor4thOrder C (
  //            2*mu + lambda, lambda, lambda,
  //            2*mu + lambda, lambda,
  //            2*mu + lambda,
  //            2*mu,
  //            2*mu,
  //            2*mu
  //            );

  tmp<volSymmTensor4thOrderField> tresult
    (
        new volSymmTensor4thOrderField
        (
            IOobject
            (
                "C",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor4thOrder
            ("C", twoMu.dimensions(), symmTensor4thOrder::zero)
        )
    );

  volSymmTensor4thOrderField& result = tresult();

  result.replace(symmTensor4thOrder::XXXX, twoMuLambda);
  result.replace(symmTensor4thOrder::XXYY, lambda);
  result.replace(symmTensor4thOrder::XXZZ, lambda);

  result.replace(symmTensor4thOrder::YYYY, twoMuLambda);
  result.replace(symmTensor4thOrder::YYZZ, lambda);

  result.replace(symmTensor4thOrder::ZZZZ, twoMuLambda);

  result.replace(symmTensor4thOrder::XYXY, twoMu);
  result.replace(symmTensor4thOrder::YZYZ, twoMu);
  result.replace(symmTensor4thOrder::ZXZX, twoMu);

  tresult().correctBoundaryConditions();

  return tresult;
}

} // End namespace Foam

// ************************************************************************* //
