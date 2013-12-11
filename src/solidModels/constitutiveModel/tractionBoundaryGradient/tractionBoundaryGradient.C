/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    tractionBoundaryGradient

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "tractionBoundaryGradient.H"
#include "fvPatch.H"
#include "Switch.H"
//#include "plasticityModel.H"
#include "constitutiveModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(tractionBoundaryGradient, 0);


// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * //

tmp<vectorField> tractionBoundaryGradient::traction
(
    const tensorField& gradField,
    const word fieldName,
    const fvPatch& patch,
    Switch orthotropic,
    word nonLinear
) const
{
    // create tmp
    tmp<vectorField> ttraction
    (
        new vectorField(gradField.size(), vector::zero)
    );
    vectorField& traction = ttraction();

    if (orthotropic)
    {
        // for now, turn off orthotropic
        FatalError
            << "tractionBoundaryGradient::traction is not written for"
            << " orthotropic yet" << nl
            << "you are probably trying to use solidContact boundaries "
            << "with orthotropic solver" << nl
            << "it should be straight-forward but I have not done it yet!"
            << exit(FatalError);
    }
    else
      {
    // lookup material properties from the solver
    const fvPatchField<scalar>& mu =
      patch.lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchField<scalar>& lambda =
      patch.lookupPatchField<volScalarField, scalar>("lambda");

    // only for nonlinear elastic properties
    // if (rheology.type() == plasticityModel::typeName)
    //   {
    //     const plasticityModel& plasticity =
    //       refCast<const plasticityModel>(rheology);

    //     mu = plasticity.newMu().boundaryField()[patch.index()];
    //     lambda = plasticity.newLambda().boundaryField()[patch.index()];
    //   }

    // required fields
    vectorField n = patch.nf();


    // calculate traction
    traction = 2*mu*(n&symm(gradField)) + lambda*tr(gradField)*n;


    //- if there is plasticity
    const constitutiveModel& rheology =
      patch.boundaryMesh().mesh().objectRegistry::
        lookupObject<constitutiveModel>("rheologyProperties");
    if (rheology.plasticityActive())
      {
        traction -=
          2*mu*(n & rheology.DEpsilonP().boundaryField()[patch.index()]);
      }

    //- if there are thermal effects
    if (patch.boundaryMesh().mesh().objectRegistry::
        foundObject<thermalModel>("thermalProperties"))
      {
        const thermalModel& thermo =
          patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<thermalModel>("thermalProperties");

        const fvPatchField<scalar>& threeKalpha =
          patch.lookupPatchField<volScalarField, scalar>
            ("((threeK*rho)*alpha)");

        if (fieldName == "DU")
          {
        const fvPatchField<scalar>& DT =
          patch.lookupPatchField<volScalarField, scalar>("DT");

        traction -=  (n*threeKalpha*DT);
          }
        else
          {
        const fvPatchField<scalar>& T =
          patch.lookupPatchField<volScalarField, scalar>("T");

        const scalarField T0 = thermo.T0()().boundaryField()[patch.index()];

        traction -=  (n*threeKalpha*(T - T0));
          }
      }

    // non linear terms
    if (nonLinear == "updatedLagrangian" || nonLinear == "totalLagrangian")
      {
        traction +=
          (n & (mu*(gradField & gradField.T())))
          + 0.5*n*lambda*(gradField && gradField);

        if (fieldName == "DU" && nonLinear == "totalLagrangian")
          {
        // incremental total Lagrangian
        const fvPatchField<tensor>& gradU =
          patch.lookupPatchField<volTensorField, tensor>("grad(U)");

        traction +=
          (n & (mu*( (gradU & gradField.T()) + (gradField & gradU.T()))))
          + 0.5*n*lambda*tr((gradU & gradField.T()) + (gradField & gradU.T()));
          }
      }


    //- add old stress for incremental approaches
    if (fieldName == "DU")
      {
        // add on old traction
        const fvPatchField<symmTensor>& sigma =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma");
        traction += (n & sigma);
      }

    //- visco-elastic
    if (rheology.viscoActive())
      {
        const fvPatchField<symmTensor>& DSigmaCorr =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("DSigmaCorr");

        traction += (n & DSigmaCorr);
      }

    //- updated Lagrangian or total Lagrangian large strain
    if (nonLinear == "updatedLagrangian" || nonLinear == "totalLagrangian")
      {
        tensorField F = I + gradField;
        if (fieldName == "DU" && nonLinear == "totalLagrangian")
          {
        const fvPatchField<tensor>& gradU =
          patch.lookupPatchField<volTensorField, tensor>("grad(U)");
        F += gradU;
          }
        tensorField Finv = inv(F);
        scalarField J = det(F);

        //- rotate second Piola Kirchhoff traction to be Cauchy traction
        traction = (traction & F)/(mag(J * Finv & n));
      }
      }

    return ttraction;
  }


// * * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * //
  tmp<vectorField> tractionBoundaryGradient::operator()
    (
     const vectorField& traction,
     const scalarField& pressure,
     const word fieldName,
     const fvPatch& patch,
     Switch orthotropic,
     word nonLinear
     ) const
  {
    // create tmp
    tmp<vectorField> tgradient(new vectorField(traction.size(), vector::zero));
    vectorField& gradient = tgradient();

    // lookup switches

    // orthotropic solvers
    if (orthotropic)
      {
    // mechanical properties
    const constitutiveModel& rheology =
      patch.boundaryMesh().mesh().objectRegistry::
        lookupObject<constitutiveModel>("rheologyProperties");
    const diagTensorField K = rheology.K()().boundaryField()[patch.index()];
    const symmTensor4thOrderField C =
        rheology.C()().boundaryField()[patch.index()];

    // required fields
    vectorField n = patch.nf();
    const diagTensorField Kinv = inv(K);
    const fvPatchField<tensor>& gradField =
      patch.lookupPatchField<volTensorField, tensor>("grad(" + fieldName + ")");


    // calculate the traction to apply
    vectorField Traction(n.size(), vector::zero);
    tensorField sigmaExp(n.size(), tensor::zero);

    //- total Lagrangian small strain
    if (fieldName == "U" && nonLinear == "off")
      {
        //- total traction
        Traction = (traction - n*pressure);

        sigmaExp = (n*(n&(C && symm(gradField)))) - (K & gradField);
      }
    //- incremental total Lagrangian small strain
    else if (fieldName == "DU" && nonLinear == "off")
      {
        const fvPatchField<symmTensor>& sigma =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma");

        //- increment of traction
        Traction = (traction - n*pressure) - (n & sigma);

        sigmaExp = (n*(n&(C && symm(gradField)))) - (K & gradField);
      }
    //- updated Lagrangian large strain
    else if (nonLinear == "updatedLagrangian")
      {
        const fvPatchField<symmTensor>& sigma =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma");

        tensorField F = I + gradField;
        tensorField Finv = inv(F);
        scalarField J = det(F);
        vectorField nCurrent = Finv & n;
        nCurrent /= mag(nCurrent);
        vectorField tractionCauchy = traction - nCurrent*pressure;

        //- increment of 2nd Piola-Kirchhoff traction
        Traction = (mag(J * Finv & n) * tractionCauchy & Finv) - (n & sigma);

        sigmaExp = (n*(n&(C && symm(gradField)))) - (K & gradField);
      }
    else
      {
        FatalError
            << "solidTractionOrtho boundary condition not suitable for "
            << " field " << fieldName << " and " << nonLinear
            << abort(FatalError);
      }

    gradient =
      n &
      (
       Kinv & ( n*(Traction) - sigmaExp )
       );
      }

    // standard isotropic solvers
    else
      {
    // lookup material properties from the solver
    const fvPatchField<scalar>& mu =
      patch.lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchField<scalar>& lambda =
      patch.lookupPatchField<volScalarField, scalar>("lambda");

    // only for nonlinear elastic properties
    // if (rheology.type() == plasticityModel::typeName)
    //   {
    //     const plasticityModel& plasticity =
    //       refCast<const plasticityModel>(rheology);

    //     mu = plasticity.newMu().boundaryField()[patch.index()];
    //     lambda = plasticity.newLambda().boundaryField()[patch.index()];
    //   }


    // required fields
    vectorField n = patch.nf();

    // gradient of the field
    const fvPatchField<tensor>& gradField =
      patch.lookupPatchField<volTensorField, tensor>("grad(" + fieldName + ")");


    //---------------------------//
    //- calculate the actual traction to apply
    //---------------------------//
    vectorField Traction(n.size(),vector::zero);

    //- total Lagrangian small strain
    if (fieldName == "U" && nonLinear == "off")
      {
        //- total traction
        Traction = (traction - n*pressure);
      }
    //- incremental total Lagrangian small strain
    else if (fieldName == "DU" && nonLinear == "off")
      {
        const fvPatchField<symmTensor>& sigma =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma");

        //- increment of traction
        Traction = (traction - n*pressure) - (n & sigma);
      }
    //- updated Lagrangian or total Lagrangian large strain
    else if (nonLinear == "updatedLagrangian" || nonLinear == "totalLagrangian")
      {
        const fvPatchField<symmTensor>& sigma =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma");

        tensorField F = I + gradField;
        if (fieldName == "DU" && nonLinear == "totalLagrangian") // for incr TL
          {
        const fvPatchField<tensor>& gradU =
          patch.lookupPatchField<volTensorField, tensor>("grad(U)");
        F += gradU;
          }
        tensorField Finv = hinv(F);
        scalarField J = det(F);
        vectorField nCurrent = Finv & n;
        nCurrent /= mag(nCurrent);
        vectorField tractionCauchy = traction - nCurrent*pressure;

        Traction = (mag(J * Finv & n) * tractionCauchy & Finv);

        if (fieldName == "DU")
          {
        //- increment of 2nd Piola-Kirchhoff traction
        Traction -= (n & sigma);
          }
      }
    else
      {
        FatalError
            << "Field " << fieldName << " and " << nonLinear
            << " nonLinear are not compatible!"
            << exit(FatalError);
      }

    //- visco-elastic
    const constitutiveModel& rheology =
      patch.boundaryMesh().mesh().objectRegistry::
        lookupObject<constitutiveModel>("rheologyProperties");
    if (rheology.viscoActive())
      {
        //Info << "visco active" << endl;
        const fvPatchField<symmTensor>& DSigmaCorr =
          patch.lookupPatchField<volSymmTensorField, symmTensor>("DSigmaCorr");

        Traction -= (n & DSigmaCorr);
      }

    //---------------------------//
    //- calculate the normal gradient based on the traction
    //---------------------------//
    gradient =
      Traction
      - (n & (mu*gradField.T() - (mu + lambda)*gradField))
      - n*lambda*tr(gradField);

    //- if there is plasticity
    if (rheology.plasticityActive())
      {
        //Info << "plasticity is active" << endl;
        gradient +=
          2*mu*(n & rheology.DEpsilonP().boundaryField()[patch.index()]);
      }

    //- if there are thermal effects
    if (patch.boundaryMesh().mesh().objectRegistry::
        foundObject<thermalModel>("thermalProperties"))
      {
        const thermalModel& thermo =
          patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<thermalModel>("thermalProperties");

        const fvPatchField<scalar>& threeKalpha =
          patch.lookupPatchField<volScalarField, scalar>
            ("((threeK*rho)*alpha)");

        if (fieldName == "DU")
          {
        const fvPatchField<scalar>& DT =
          patch.lookupPatchField<volScalarField, scalar>("DT");

        gradient +=  (n*threeKalpha*DT);
          }
        else
          {
        const fvPatchField<scalar>& T =
          patch.lookupPatchField<volScalarField, scalar>("T");

        const scalarField T0 = thermo.T0()().boundaryField()[patch.index()];

        gradient +=  (n*threeKalpha*(T - T0));
          }
      }

    //- higher order non-linear terms
    if (nonLinear == "updatedLagrangian" || nonLinear == "totalLagrangian")
      {
          // no extra relaxation
          gradient -=
              (n & (mu*(gradField & gradField.T())))
              // + 0.5*n*lambda*(gradField && gradField);
              + 0.5*n*lambda*tr(gradField & gradField.T());
          //- tensorial identity
          //- tr(gradField & gradField.T())*I == (gradField && gradField)*I

          if (fieldName == "DU" && nonLinear == "totalLagrangian")
          {
              // gradU is const in a time step
              const fvPatchField<tensor>& gradU =
                  patch.lookupPatchField<volTensorField, tensor>("grad(U)");

              gradient -=
                  (n &
                   (mu*( (gradField & gradU.T())
                         + (gradU & gradField.T()) ))
                      )
                  + 0.5*n*lambda*tr( (gradField & gradU.T())
                                     + (gradU & gradField.T()) );
          }
      }

    gradient /= (2.0*mu + lambda);
      }

    return tgradient;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
