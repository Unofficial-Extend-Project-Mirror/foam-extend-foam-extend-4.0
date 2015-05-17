    label cohesivePatchID = -1;

    solidCohesiveFvPatchVectorField* cohesivePatchUPtr = NULL;
    solidCohesiveFixedModeMixFvPatchVectorField* cohesivePatchUFixedModePtr = NULL;

    forAll (U.boundaryField(), patchI)
    {
        if (isA<solidCohesiveFvPatchVectorField>(U.boundaryField()[patchI]))
        {
            cohesivePatchID = patchI;
            cohesivePatchUPtr =
                &refCast<solidCohesiveFvPatchVectorField>
                (
                     U.boundaryField()[cohesivePatchID]
                );
            break;
        }
        else if (isA<solidCohesiveFixedModeMixFvPatchVectorField>(U.boundaryField()[patchI]))
        {
            cohesivePatchID = patchI;
            cohesivePatchUFixedModePtr =
                &refCast<solidCohesiveFixedModeMixFvPatchVectorField>
                (
                    U.boundaryField()[cohesivePatchID]
                );
            break;
        }
    }

    if(cohesivePatchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Can't find cohesiveLawFvPatch" << nl
            << "One of the boundary patches in " << U.name() << ".boundaryField() "
            << "should be of type " << solidCohesiveFvPatchVectorField::typeName
            << "or " << solidCohesiveFixedModeMixFvPatchVectorField::typeName
            << abort(FatalError);
    }

    // solidCohesiveFvPatchVectorField& cohesivePatchU =
    //     refCast<solidCohesiveFvPatchVectorField>
    //     (
    //         U.boundaryField()[cohesivePatchID]
    //     );

    // philipc: I have moved cohesive stuff to constitutiveModel
    // cohesiveZone is an index field
    // which allows the user to limit the crack to certain areas at runtime
    // 1 for faces within cohesiveZone
    // 0 for faces outside cohesiveZone
    surfaceScalarField cohesiveZone
    (
        IOobject
        (
            "cohesiveZone",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
     );

    // limit crack to specified boxes
    {
        const dictionary& stressControl =
            mesh.solutionDict().subDict("solidMechanics");

        List<boundBox> userBoxes(stressControl.lookup("crackLimitingBoxes"));
        const surfaceVectorField& Cf = mesh.Cf();
        forAll(cohesiveZone.internalField(), faceI)
        {
            bool faceInsideBox = false;

            forAll(userBoxes, boxi)
            {
                if(userBoxes[boxi].contains(Cf.internalField()[faceI])) faceInsideBox = true;
            }

            if(faceInsideBox)
            {
                cohesiveZone.internalField()[faceI] = 1.0;
            }
        }

        forAll(cohesiveZone.boundaryField(), patchI)
        {
            // cracks may go along proc boundaries
            if(mesh.boundaryMesh()[patchI].type() == processorPolyPatch::typeName)
            {
                forAll(cohesiveZone.boundaryField()[patchI], faceI)
                {
                    bool faceInsideBox = false;

                    forAll(userBoxes, boxi)
                    {
                        if(userBoxes[boxi].contains(Cf.boundaryField()[patchI][faceI])) faceInsideBox = true;
                    }

                    if(faceInsideBox)
                    {
                        cohesiveZone.boundaryField()[patchI][faceI] = 1.0;
                    }
                }
            }
        }

        Info << "\nThere are " << gSum(cohesiveZone.internalField()) << " potential internal crack faces" << nl << endl;
        Info << "\nThere are " << gSum(cohesiveZone.boundaryField())/2 << " potential coupled boundary crack faces" << nl << endl;

        // write field for visualisation
        volScalarField cohesiveZoneVol
        (
            IOobject
            (
                "cohesiveZoneVol",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        forAll(cohesiveZone.internalField(), facei)
        {
            if(cohesiveZone.internalField()[facei])
            {
                cohesiveZoneVol.internalField()[mesh.owner()[facei]] = 1.0;
                cohesiveZoneVol.internalField()[mesh.neighbour()[facei]] = 1.0;
            }
        }

        forAll(cohesiveZone.boundaryField(), patchi)
        {
            forAll(cohesiveZone.boundaryField()[patchi], facei)
            {
                if(cohesiveZone.boundaryField()[patchi][facei] > 0.0)
                {
                    cohesiveZoneVol.boundaryField()[patchi][facei] = 1.0;
                }
            }
        }
        Info << "Writing cohesiveZone field" << endl;
        cohesiveZoneVol.write();
    }


    Switch initialiseSolution(false);

    if
    (
        mesh.solutionDict().subDict("solidMechanics")
       .found("initialiseSolution")
    )
    {
        initialiseSolution =
            Switch
            (
                mesh.solutionDict().subDict("solidMechanics").lookup
                (
                    "initialiseSolution"
                )
            );
    }


    Switch breakOnlyOneFacePerTopologyChange(true);

    if
    (
        mesh.solutionDict().subDict("solidMechanics")
       .found("breakOnlyOneFacePerTopologyChange")
    )
    {
        breakOnlyOneFacePerTopologyChange =
            Switch
            (
                mesh.solutionDict().subDict("solidMechanics").lookup
                (
                    "breakOnlyOneFacePerTopologyChange"
                )
            );
    }


    Switch crackPropagationFromSpecifiedPatches
    (
        mesh.solutionDict().subDict("solidMechanics").lookup
        (
            "crackPropagationFromSpecifiedPatches"
        )
    );

    wordList crackPropagationPatchNames
    (
        mesh.solutionDict().subDict("solidMechanics").lookup
        (
            "crackPropagationPatches"
        )
    );

    labelList crackPropagationPatches(crackPropagationPatchNames.size(), -1);

    forAll(crackPropagationPatchNames, patchI)
    {
        crackPropagationPatches[patchI] =
            mesh.boundaryMesh().findPatchID
            (
                crackPropagationPatchNames[patchI]
            );

        if(crackPropagationPatches[patchI] == -1)
        {
            FatalErrorIn(args.executable())
                << "Can't find " << crackPropagationPatchNames[patchI]
                << " patch" << abort(FatalError);
        }
    }

    // Internal faces next to selected crack propagation patches
    labelList crackPropagationPatchesInternalFaces;

#   include "updateCrackPropagationPatchesInternalFaces.H"

