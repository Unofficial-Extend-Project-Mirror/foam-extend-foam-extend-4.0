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


// * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::tractionBoundaryGradient::traction
(
    const tensorField& gradField,
    const word& workingFieldName,
    const word& integralFieldName,
    const fvPatch& patch,
    const bool orthotropic,
    const nonLinearGeometry::nonLinearType& nonLinear,
    const bool incremental
)
{
    // Create result
    tmp<vectorField> ttraction
    (
        new vectorField(gradField.size(), vector::zero)
    );
    vectorField& traction = ttraction();

    // Orthotropic material
    if (orthotropic)
    {
        // For now, turn off orthotropic
        FatalErrorIn
        (
            "tmp<vectorField> tractionBoundaryGradient::traction\n"
            "(\n"
            "    const tensorField& gradField,\n"
            "    const word& workingFieldName,\n"
            "    const word& integralFieldName,\n"
            "    const fvPatch& patch,\n"
            "    const bool orthotropic,\n"
            "    const nonLinearGeometry::nonLinearType& nonLinear\n"
            ") const"
        )   << "tractionBoundaryGradient::traction is not written for"
            << " orthotropic materials yet" << nl
            << "you are probably trying to use solidContact boundaries "
            << "with the orthotropic solver" << nl
            << exit(FatalError);
    }
    else
    {
        // Lookup material properties from the solver
        const fvPatchScalarField& mu =
            patch.lookupPatchField<volScalarField, scalar>("mu");

        const fvPatchScalarField& lambda =
            patch.lookupPatchField<volScalarField, scalar>("lambda");

        // only for nonlinear elastic properties
        // if (rheology.type() == plasticityModel::typeName)
        // {
        //     const plasticityModel& plasticity =
        //       refCast<const plasticityModel>(rheology);

        //     mu = plasticity.newMu().boundaryField()[patch.index()];
        //     lambda = plasticity.newLambda().boundaryField()[patch.index()];
        // }

        // Get patch normal
        vectorField n = patch.nf();

        // Calculate traction
        traction = 2*mu*(n & symm(gradField)) + lambda*tr(gradField)*n;

        // Plasticity effects
        const constitutiveModel& rheology =
            patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<constitutiveModel>("rheologyProperties");

        if (rheology.plasticityActive())
        {
            traction -=
                2*mu*(n & rheology.DEpsilonP().boundaryField()[patch.index()]);
        }

        // Thermal effects
        if
        (
            patch.boundaryMesh().mesh().objectRegistry::
            foundObject<thermalModel>("thermalProperties")
        )
        {
            const thermalModel& thermo =
                patch.boundaryMesh().mesh().objectRegistry::
                lookupObject<thermalModel>("thermalProperties");

            const fvPatchScalarField& threeKalpha =
                patch.lookupPatchField<volScalarField, scalar>
                ("((threeK*rho)*alpha)");

            // Incremental thermal contribution
            if (incremental)
            {
                const fvPatchScalarField& DT =
                    patch.lookupPatchField<volScalarField, scalar>("DT");

                traction -= (n*threeKalpha*DT);
            }
            else
            {
                const fvPatchScalarField& T =
                    patch.lookupPatchField<volScalarField, scalar>("T");

                const scalarField T0 =
                    thermo.T0()().boundaryField()[patch.index()];

                traction -= (n*threeKalpha*(T - T0));
            }
        }

        // Non-linear terms
        if
        (
            nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            traction +=
                (n & (mu*(gradField & gradField.T())))
              + 0.5*n*lambda*(gradField && gradField);

            if
            (
                incremental
             && nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
            )
            {
                // Incremental total Lagrangian

                const fvPatchTensorField& gradU =
                    patch.lookupPatchField<volTensorField, tensor>
                    (
                        "grad(" + integralFieldName + ")"
                    );

                traction +=
                    (
                        n
                      & (
                          mu*
                          (
                              (gradU & gradField.T())
                            + (gradField & gradU.T())
                            )
                        )
                    )
                  + 0.5*n*lambda*tr
                    (
                        (gradU & gradField.T())
                      + (gradField & gradU.T())
                    );
            }
        }


        // Add old stress for incremental approaches
        if (incremental)
        {
            // add on old traction
            const fvPatchSymmTensorField& sigma =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigma"
                );

            traction += (n & sigma);
        }

        // Visco-elastic effects
        if (rheology.viscoActive())
        {
            const fvPatchSymmTensorField& DSigmaCorr =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "DSigmaCorr"
                );

            traction += (n & DSigmaCorr);
        }

        // Updated Lagrangian or total Lagrangian large strain
        if
        (
            nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            tensorField F = I + gradField;

            if
            (
                incremental
             && nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
            )
            {
                // Incremental total Lagrangian

                const fvPatchTensorField& gradU =
                    patch.lookupPatchField<volTensorField, tensor>
                    (
                        "grad(" + integralFieldName + ")"
                    );

                F += gradU;
            }

            tensorField Finv = inv(F);
            scalarField J = det(F);

            // Rotate second Piola Kirchhoff traction to be Cauchy traction
            traction = (traction & F)/(mag(J * Finv & n));
        }
    }

    return ttraction;
}


// * * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::tractionBoundaryGradient::snGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const word& workingFieldName,
    const word& integralFieldName,
    const fvPatch& patch,
    const bool orthotropic,
    const nonLinearGeometry::nonLinearType& nonLinear,
    const bool incremental
)
{
    // Create result
    tmp<vectorField> tgradient(new vectorField(traction.size(), vector::zero));
    vectorField& gradient = tgradient();

    // Orthotropic material
    if (orthotropic)
    {
        // Get mechanical properties
        const constitutiveModel& rheology =
            patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<constitutiveModel>("rheologyProperties");

        const diagTensorField K =
            rheology.K()().boundaryField()[patch.index()];

        const symmTensor4thOrderField C =
            rheology.C()().boundaryField()[patch.index()];

        // Required fields
        vectorField n = patch.nf();
        const diagTensorField Kinv = inv(K);
        const fvPatchTensorField& gradField =
            patch.lookupPatchField<volTensorField, tensor>
            (
                "grad(" + workingFieldName + ")"
            );

        // Calculate the traction to apply
        vectorField Traction(n.size(), vector::zero);
        tensorField sigmaExp(n.size(), tensor::zero);

        // Total Lagrangian, small strain
        if
        (
            !incremental
         && nonLinear == nonLinearGeometry::OFF
        )
        {
            // Use total traction
            Traction = (traction - n*pressure);

            sigmaExp = (n*(n&(C && symm(gradField)))) - (K & gradField);
        }
        // Incremental total Lagrangian small strain
        else if
        (
            incremental
         && nonLinear == nonLinearGeometry::OFF
        )
        {
            const fvPatchSymmTensorField& sigma =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigma"
                );

            // Increment of traction
            Traction = (traction - n*pressure) - (n & sigma);

            sigmaExp = (n*(n&(C && symm(gradField)))) - (K & gradField);
        }
        // Updated Lagrangian large strain
        else if
        (
            nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
        )
        {
            const fvPatchSymmTensorField& sigma =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigma"
                );

            tensorField F = I + gradField;
            tensorField Finv = inv(F);
            scalarField J = det(F);
            vectorField nCurrent = Finv & n;
            nCurrent /= mag(nCurrent) + SMALL;
            vectorField tractionCauchy = traction - nCurrent*pressure;

            // Increment of 2nd Piola-Kirchhoff traction
            Traction = mag(J*(Finv & n))*(tractionCauchy & Finv) - (n & sigma);

            sigmaExp = n*(n &(C && symm(gradField))) - (K & gradField);
        }
        else
        {
            FatalErrorIn
            (
                "tmp<vectorField> tractionBoundaryGradient::snGrad\n"
                "(\n"
                "    const vectorField& traction,\n"
                "    const scalarField& pressure,\n"
                "    const word& workingFieldName,\n"
                "    const word& integralFieldName,\n"
                "    const fvPatch& patch,\n"
                "    const bool orthotropic,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinear,\n"
                "    const bool incremental\n"
                ") const"
            )   << "solidTractionOrtho boundary condition not suitable for "
                << " field " << workingFieldName << " and "
                << nonLinearGeometry::nonLinearNames_[nonLinear]
                << abort(FatalError);
        }

        gradient = n & (Kinv & ( n*(Traction) - sigmaExp ));
    }
    else
    {
        // Standard isotropic solvers

        // Lookup material properties from the solver
        const fvPatchScalarField& mu =
            patch.lookupPatchField<volScalarField, scalar>("mu");

        const fvPatchScalarField& lambda =
            patch.lookupPatchField<volScalarField, scalar>("lambda");

        // only for nonlinear elastic properties
        // if (rheology.type() == plasticityModel::typeName)
        // {
        //     const plasticityModel& plasticity =
        //       refCast<const plasticityModel>(rheology);

        //     mu = plasticity.newMu().boundaryField()[patch.index()];
        //     lambda = plasticity.newLambda().boundaryField()[patch.index()];
        // }

        vectorField n = patch.nf();

        // gradient of the field
        const fvPatchTensorField& gradField =
            patch.lookupPatchField<volTensorField, tensor>
            (
                "grad(" + workingFieldName + ")"
            );


        vectorField Traction(n.size(), vector::zero);

        // Total Lagrangian, small strain
        if
        (
            !incremental
         && nonLinear == nonLinearGeometry::OFF
        )
        {
            // Use total traction
            Traction = (traction - n*pressure);
        }
        // Incremental total Lagrangian small strain
        else if
        (
            incremental
         && nonLinear == nonLinearGeometry::OFF
        )
        {
            const fvPatchSymmTensorField& sigma =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigma"
                );

            // Increment of traction
            Traction = (traction - n*pressure) - (n & sigma);
        }
        else if
        (
            nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            const fvPatchSymmTensorField& sigma =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigma"
                );

            tensorField F = I + gradField;

            if
            (
                incremental
             && nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
            )
            {
                // Incremental total Lagrangian

                const fvPatchTensorField& gradU =
                    patch.lookupPatchField<volTensorField, tensor>
                    (
                        "grad(" + integralFieldName + ")"
                    );

                F += gradU;
            }

            tensorField Finv = hinv(F);
            scalarField J = det(F);
            vectorField nCurrent = Finv & n;
            nCurrent /= mag(nCurrent) + SMALL;
            vectorField tractionCauchy = traction - nCurrent*pressure;

            Traction = mag(J*(Finv & n))*(tractionCauchy & Finv);

            if (incremental)
            {
                // Increment of 2nd Piola-Kirchhoff traction
                Traction -= (n & sigma);
            }
        }
        else
        {
            FatalErrorIn
            (
                "tmp<vectorField> tractionBoundaryGradient::snGrad\n"
                "(\n"
                "    const vectorField& traction,\n"
                "    const scalarField& pressure,\n"
                "    const word& workingFieldName,\n"
                "    const word& integralFieldName,\n"
                "    const fvPatch& patch,\n"
                "    const bool orthotropic,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinear,\n"
                "    const bool incremental\n"
                ") const"
            )   << " field " << workingFieldName << " and "
                << nonLinearGeometry::nonLinearNames_[nonLinear]
                << " nonLinear are not compatible!"
                << abort(FatalError);
      }

        //- visco-elastic
        const constitutiveModel& rheology =
            patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<constitutiveModel>("rheologyProperties");

        if (rheology.viscoActive())
        {
            const fvPatchSymmTensorField& DSigmaCorr =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "DSigmaCorr"
                );

            Traction -= (n & DSigmaCorr);
        }

        // Calculate the normal gradient based on the traction

        gradient =
            Traction
            - (n & (mu*gradField.T() - (mu + lambda)*gradField))
            - n*lambda*tr(gradField);

        //- Plasticity contribution
        if (rheology.plasticityActive())
        {
            gradient +=
                2*mu*(n & rheology.DEpsilonP().boundaryField()[patch.index()]);
        }

        // Thermal effects
        if
        (
            patch.boundaryMesh().mesh().objectRegistry::
            foundObject<thermalModel>("thermalProperties")
        )
        {
            const thermalModel& thermo =
                patch.boundaryMesh().mesh().objectRegistry::
                lookupObject<thermalModel>("thermalProperties");

            const fvPatchScalarField& threeKalpha =
                patch.lookupPatchField<volScalarField, scalar>
                (
                    "((threeK*rho)*alpha)"
                );

                if (!incremental)
                {
                    const fvPatchScalarField& DT =
                        patch.lookupPatchField<volScalarField, scalar>("DT");

                    gradient += n*threeKalpha*DT;
                }
                else
                {
                    const fvPatchScalarField& T =
                        patch.lookupPatchField<volScalarField, scalar>("T");

                    const scalarField T0 =
                        thermo.T0()().boundaryField()[patch.index()];

                    gradient += n*threeKalpha*(T - T0);
                }
        }

        // Higher order non-linear terms
        if
        (
            nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            // no extra relaxation
            gradient -=
                (n & (mu*(gradField & gradField.T())))
                // + 0.5*n*lambda*(gradField && gradField);
              + 0.5*n*lambda*tr(gradField & gradField.T());

            if
            (
                incremental
             && nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
            )
            {
                // gradU is const in a time step
                const fvPatchTensorField& gradU =
                    patch.lookupPatchField<volTensorField, tensor>("grad(U)");

                gradient -=
                    (
                        n &
                        (
                            mu*
                            (
                                (gradField & gradU.T())
                              + (gradU & gradField.T())
                            )
                        )
                    )
                  + 0.5*n*lambda*tr
                    (
                        (gradField & gradU.T())
                      + (gradU & gradField.T())
                    );
            }
        }

        gradient /= (2.0*mu + lambda);
    }

    return tgradient;
}


// ************************************************************************* //
