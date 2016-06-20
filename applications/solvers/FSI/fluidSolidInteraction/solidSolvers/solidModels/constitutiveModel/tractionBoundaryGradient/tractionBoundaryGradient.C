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
#include "constitutiveModel.H"
#include "thermalModel.H"


// * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::tractionBoundaryGradient::traction
(
    const tensorField& gradField,
    const word& workingFieldName,
    const word& integralFieldName,
    const fvPatch& patch,
    const bool incremental
)
{
    // Lookup nonLinear type from solver
    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            patch.boundaryMesh().mesh().solutionDict().subDict
            (
                "solidMechanics"
            ).lookup("nonLinear")
        );

//     Info << nonLinearGeometry::nonLinearNames_ << endl;

    // Lookup if solver uses an orthotropic approach
    const Switch orthotropic =
        patch.boundaryMesh().mesh().solutionDict().subDict
        (
            "solidMechanics"
        ).lookupOrDefault<Switch>("orthotropic", false);

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
        // if (mechanical.type() == plasticityModel::typeName)
        // {
        //     const plasticityModel& plasticity =
        //       refCast<const plasticityModel>(mechanical);

        //     mu = plasticity.newMu().boundaryField()[patch.index()];
        //     lambda = plasticity.newLambda().boundaryField()[patch.index()];
        // }

        // Get patch normal
        vectorField n = patch.nf();

        // Calculate traction
        traction = 2*mu*(n & symm(gradField)) + lambda*tr(gradField)*n;

        const fvMesh& mesh = patch.boundaryMesh().mesh();

        // Plasticity effects
        const constitutiveModel& mechanical =
            mesh.lookupObject<constitutiveModel>("rheologyProperties");

        if (mesh.foundObject<volSymmTensorField>("DEpsilonP"))
        {
            const fvPatchSymmTensorField& DEpsilonP =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "DEpsilonP"
                );

            traction -= 2*mu*(n & DEpsilonP) + lambda*n*tr(DEpsilonP);
        }

        // Thermal effects
        if
        (
            patch.boundaryMesh().mesh().objectRegistry::
            foundObject<thermalModel>("thermalProperties")
        )
        {
            // const thermalModel& thermo =
            //    patch.boundaryMesh().mesh().objectRegistry::
            //    lookupObject<thermalModel>("thermalProperties");

            if
            (
                patch.boundaryMesh().mesh().objectRegistry::
                foundObject<volScalarField>("(threeK*alpha)")
            )
            {
                const fvPatchScalarField& threeKalpha =
                    patch.lookupPatchField<volScalarField, scalar>
                    (
                        "(threeK*alpha)"
                    );

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

                //const scalarField T0 =
                //    thermo.T0()().boundaryField()[patch.index()];

                // traction -= (n*threeKalpha*(T - T0));
                    traction -= (n*threeKalpha*T);
                }
            }
            else
            {
                const fvPatchScalarField& threeK =
                    patch.lookupPatchField<volScalarField, scalar>
                    (
                        "threeK"
                    );

                const fvPatchScalarField& alpha =
                    patch.lookupPatchField<volScalarField, scalar>
                    (
                        "alpha"
                    );

                const fvPatchScalarField& DT =
                    patch.lookupPatchField<volScalarField, scalar>("DT");
                traction -= (n*threeK*alpha*DT);

                // // Incremental thermal contribution
                // if (incremental)
                // {

                //     traction -= (n*threeK*alpha*DT);
                // }
                // else
                // {
                //     const fvPatchScalarField& T =
                //         patch.lookupPatchField<volScalarField, scalar>("T");

                // const scalarField T0 =
                //     thermo.T0()().boundaryField()[patch.index()];
                // traction -= (n*threeK*alpha*(T - T0));

                //     // traction -= (n*threeK*alpha*T);
                // }
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
        if (mechanical.viscoActive())
        {
            if (workingFieldName == "DU" || workingFieldName == "DD")
            {
                //Info << "visco active" << endl;
                const fvPatchField<symmTensor>& DSigmaCorr =
                    patch.lookupPatchField<volSymmTensorField, symmTensor>
                    (
                        "DSigmaCorr"
                    );

                traction -= n & DSigmaCorr;
            }
            else if (workingFieldName == "U" || workingFieldName == "D")
            {
                //Info << "visco active" << endl;
                const fvsPatchField<symmTensor>& sigmaCorrf =
                    patch.lookupPatchField
                    <
                        surfaceSymmTensorField,
                        symmTensor
                    >
                    (
                        "sigmaCorrf"
                    );

                traction -= n & sigmaCorrf;
            }
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
        else if (nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF)
        {
            if (!incremental)
            {
                FatalError
                    << "tractionBoundaryGradient.traction() "
                    << "field name must be DU with updated Lagrangian Kirchhoff"
                    << abort(FatalError);
            }

            if (mechanical.viscoActive())
            {
                FatalError
                    << "tractionBoundaryGradient.traction() "
                    << "nonLinear updatedLagrangianKirchhoff with viscous"
                    << " stresses have not been implemented yet"
                    << abort(FatalError);
            }

            // Lookup Tau from the solver as it is up-to-date
            const fvPatchField<symmTensor>& tau =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "tauKirchhoff"
                );
            // Relative inverse deformation gradient
            const fvPatchField<tensor>& Finv =
                patch.lookupPatchField<volTensorField, tensor>
                (
                    "relFinv"
                );
            // Surface fields
            // const fvsPatchField<symmTensor>& tau =
            //     patch.lookupPatchField<surfaceSymmTensorField, symmTensor>
            //     (
            //         "tauKirchhofff"
            //     );
            // const fvsPatchField<tensor>& Finv =
            //     patch.lookupPatchField<surfaceTensorField, symmTensor>
            //     (
            //         "relFinvf"
            //     );

            // Calculate traction with deformed normal
            // Note that J drops out as sigmaCauchy = tau/J
            traction = (Finv.T() & n) & tau;
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
    const bool incremental
)
{
    // Lookup nonLinear type from solver
    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            patch.boundaryMesh().mesh().solutionDict().subDict
            (
                "solidMechanics"
            ).lookup("nonLinear")
        );

    // Lookup if solver uses an orthotropic approach
    const Switch orthotropic =
        patch.boundaryMesh().mesh().solutionDict().subDict
        (
            "solidMechanics"
        ).lookupOrDefault<Switch>("orthotropic", false);

    // Create result
    tmp<vectorField> tgradient(new vectorField(traction.size(), vector::zero));
    vectorField& gradient = tgradient();

    // Orthotropic material
    if (orthotropic)
    {
        // Get mechanical properties
        const constitutiveModel& mechanical =
            patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<constitutiveModel>("rheologyProperties");

        const diagTensorField K =
            mechanical.K()().boundaryField()[patch.index()];

        const symmTensor4thOrderField C =
            mechanical.C()().boundaryField()[patch.index()];

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

            sigmaExp = (n*(n & (C && symm(gradField)))) - (K & gradField);
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

            sigmaExp = (n*(n & (C && symm(gradField)))) - (K & gradField);
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

            sigmaExp = n*(n & (C && symm(gradField))) - (K & gradField);
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
                "    const bool incremental\n"
                ") const"
            )   << "solidTractionOrtho boundary condition not suitable for "
                << " field " << workingFieldName << " and "
                << nonLinearGeometry::nonLinearNames_[nonLinear]
                << abort(FatalError);
        }
    }
    else
    {
        // Standard isotropic solvers

        // Lookup material properties from the solver
        const fvPatchScalarField& mu =
            patch.lookupPatchField<volScalarField, scalar>("mu");

        const fvPatchScalarField& lambda =
            patch.lookupPatchField<volScalarField, scalar>("lambda");

        // Patch unit normals
        vectorField n = patch.nf();

        // Lookup gradient of the working field
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
            // sigma: sigma.oldTime
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
        else if (nonLinear == nonLinearGeometry::DEFORMED_LAGRANGIAN)
        {
            if (!incremental)
            {
                FatalError
                    << "tractionBoundaryGradient: Field " << workingFieldName
                    << " and deformed Lagranngian  are not compatible!"
                    << abort(FatalError);
            }

            Traction = traction - n*pressure;

            const fvPatchField<symmTensor>& oldSigmaCauchy =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigmaCauchy"
                );
            Traction -= n & oldSigmaCauchy;

            // const fvsPatchField<symmTensor>& rotatedSigmaOld =
            //     patch.lookupPatchField<
            //         surfaceSymmTensorField, symmTensor
            //         >("rotatedSigmaOldf");
            // Traction -= n & rotatedSigmaOld;
        }

        // Visco-elastic
        const constitutiveModel& mechanical =
            patch.boundaryMesh().mesh().objectRegistry::
            lookupObject<constitutiveModel>("rheologyProperties");

        if (mechanical.viscoActive())
        {
            if (workingFieldName == "DU" || workingFieldName == "DD")
            {
                //Info << "visco active" << endl;
                const fvPatchField<symmTensor>& DSigmaCorr =
                    patch.lookupPatchField<volSymmTensorField, symmTensor>
                    (
                        "DSigmaCorr"
                    );

                Traction -= n & DSigmaCorr;
            }
            else if (workingFieldName == "U" || workingFieldName == "D")
            {
                //Info << "visco active" << endl;
                const fvsPatchField<symmTensor>& sigmaCorrf =
                    patch.lookupPatchField
                    <
                        surfaceSymmTensorField,
                    symmTensor
                    >
                    (
                        "sigmaCorrf"
                    );

                Traction -= n & sigmaCorrf;
            }
        }

        // Calculate the normal gradient based on the traction
        gradient =
            Traction
            - (n & (mu*gradField.T() - (mu + lambda)*gradField))
            - n*lambda*tr(gradField);

        const fvMesh& mesh = patch.boundaryMesh().mesh();

        // Plasticity contribution
        if (mesh.foundObject<volSymmTensorField>("DEpsilonP"))
        {
            const fvPatchSymmTensorField& DEpsilonP =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "DEpsilonP"
                );

            gradient += 2*mu*(n & DEpsilonP) + lambda*n*tr(DEpsilonP);
        }

        // Thermal effects
        if
        (
            patch.boundaryMesh().mesh().objectRegistry::
            foundObject<thermalModel>("thermalProperties")
        )
        {
            // const thermalModel& thermo =
            //    patch.boundaryMesh().mesh().objectRegistry::
            //    lookupObject<thermalModel>("thermalProperties");

            if
            (
                patch.boundaryMesh().mesh().objectRegistry::
                foundObject<volScalarField>("(threeK*alpha)")
            )
            {
                const fvPatchScalarField& threeKalpha =
                    patch.lookupPatchField<volScalarField, scalar>
                    (
                        "(threeK*alpha)"
                    );

                // if (!incremental)
                if (incremental) // bugfix philipc Nov-14
                {
                    const fvPatchScalarField& DT =
                        patch.lookupPatchField<volScalarField, scalar>
                        (
                           "DT"
                        );

                    gradient += n*threeKalpha*DT;
                }
                else
                {
                    const fvPatchScalarField& T =
                        patch.lookupPatchField<volScalarField, scalar>
                        (
                            "T"
                        );

                    const fvPatchScalarField& T0 =
                        patch.lookupPatchField<volScalarField, scalar>
                        (
                            "T0"
                        );

                //const scalarField T0 =
                //    thermo.T0()().boundaryField()[patch.index()];
                //gradient += n*threeKalpha*T;

                    gradient += n*threeKalpha*(T - T0);
                }
            }
            else
            {
                const fvPatchScalarField& threeK =
                    patch.lookupPatchField<volScalarField, scalar>
                    (
                        "threeK"
                    );

                const fvPatchScalarField& alpha =
                    patch.lookupPatchField<volScalarField, scalar>
                    (
                        "alpha"
                    );

                const fvPatchScalarField& DT =
                    patch.lookupPatchField<volScalarField, scalar>("DT");
                gradient += (n*threeK*alpha*DT);
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
        else if (nonLinear == nonLinearGeometry::DEFORMED_LAGRANGIAN)
        {
            if (!incremental)
            {
                FatalError
                    << "tractionBoundaryGradient: Field " << workingFieldName
                    << " and deformed Lagranngian  are not compatible!"
                    << abort(FatalError);
            }

            const fvPatchField<tensor>& skewGradDU =
                patch.lookupPatchField<volTensorField, tensor>
                (
                    "skew(gradDU)"
                );
            const fvPatchField<symmTensor>& sigmaCauchy =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigmaCauchy"
                );

            gradient +=
                (n & (mu*(gradField & gradField.T())))
              + 0.5*n*lambda*tr(gradField & gradField.T())
              - (
                    n &
                    (
                        (skewGradDU & sigmaCauchy) - (sigmaCauchy & skewGradDU)
                    )
                );
        }
        else if (nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF)
        {
            // Relative deformation gradient inverse
            const fvPatchField<tensor>& Finv =
                 patch.lookupPatchField<volTensorField, tensor>("relFinv");
            // tau is up-to-date in the solver
            const fvPatchField<symmTensor>& tau =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "tauKirchhoff"
                );
            // Surface fields
            // const fvsPatchField<tensor>& Finv =
            //     patch.lookupPatchField<surfaceTensorField, tensor>
            //     (
            //         "relFinvf"
            //     );
            // const fvsPatchField<symmTensor>& tau =
            //     patch.lookupPatchField<surfaceSymmTensorField, symmTensor>
            //     (
            //         "tauKirchhofff"
            //     );
            // Total Jacobian relates tau and sigmaCauchy
            // const fvPatchField<scalar>& J =
            //      patch.lookupPatchField<volScalarField, scalar>("relJ");
            const fvPatchField<scalar>& J =
                patch.lookupPatchField<volScalarField, scalar>("J");

            vectorField nCurrent = J*Finv.T() & n;
            // nCurrent should be unit vectors but we will force normalisation
            // to remove any errors
            nCurrent /= mag(nCurrent);
            vectorField tractionCauchy = (traction - nCurrent*pressure);

            // We overwrite the gradient
            gradient =
                tractionCauchy
                - (nCurrent & tau/J)
                + (2.0*mu + lambda)*(n & gradField);
        }
        else if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN_KIRCHHOFF)
        {
            // Relative deformation gradient inverse
            const fvPatchField<tensor>& Finv =
                patch.lookupPatchField<volTensorField, tensor>("Finv");
            // tau is up-to-date in the solver
            const fvPatchField<symmTensor>& tau =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "tauKirchhoff"
                );
            // Total Jacobian relates tau and sigmaCauchy
            const fvPatchField<scalar>& J =
                patch.lookupPatchField<volScalarField, scalar>("J");
            // const fvPatchField<scalar>& J =
            //      patch.lookupPatchField<volScalarField, scalar>("relJ");

            vectorField nCurrent = J*Finv.T() & n;
            // nCurrent should be unit vectors but we will force normalisation
            // to remove any errors
            nCurrent /= mag(nCurrent);
            vectorField tractionCauchy = (traction - nCurrent*pressure);

            gradient =
                tractionCauchy
                - (nCurrent & tau/J)
                + (2.0*mu + lambda)*(n & gradField);
        }
        // else if (nonLinear == nonLinearGeometry::ANISOTROPIC_PLASTICITY)
        // {
        //     // Over-write previous gradient calculations
        //     // change sigma to whatever the stress is called in the solver

        //     // Lookup sigma from the solver as it is up-to-date inside the
        //     // momentum loop
        //     const fvPatchField<symmTensor>& sigmaCauchy =
        //         patch.lookupPatchField<volSymmTensorField, symmTensor>
        //         (
        //             "sigmaCauchy"
        //         );

        //     gradient =
        //         traction - (n & sigmaCauchy)
        //         + (2.0*mu + lambda)*(n & gradField);
        // }

        gradient /= (2.0*mu + lambda);
    }

    return tgradient;
}


// ************************************************************************* //
