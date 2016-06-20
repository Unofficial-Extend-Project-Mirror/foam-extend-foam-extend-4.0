{
    // DEpsilonPEq
    forAll(DEpsilonPEqf_.boundaryField(), patchI)
    {
        if
        (
            isA<processorFvsPatchScalarField>
            (
                DEpsilonPEqf_.boundaryField()[patchI]
            )
        )
        {
            if (Pstream::parRun())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                const scalarField& f =
                    DEpsilonPEqf_.boundaryField()[patchI];

                OPstream::write
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(f.begin()),
                    f.byteSize()
                );
            }
        }
    }

    forAll(DEpsilonPEqf_.boundaryField(), patchI)
    {
        if
        (
            isA<processorFvsPatchScalarField>
            (
                DEpsilonPEqf_.boundaryField()[patchI]
            )
        )
        {
            if (Pstream::parRun())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                ngbProcDEpsilonPEq.set
                (
                    patchI,
                    new scalarField
                    (
                        DEpsilonPEqf_.boundaryField()[patchI].size(),
                        0
                    )
                );

                IPstream::read
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>
                    (
                        ngbProcDEpsilonPEq[patchI].begin()
                    ),
                    ngbProcDEpsilonPEq[patchI].byteSize()
                );
            }
        }
    }

    // DEpsilonP
    forAll(DEpsilonPf_.boundaryField(), patchI)
    {
        if
        (
            isA<processorFvsPatchSymmTensorField>
            (
                DEpsilonPf_.boundaryField()[patchI]
            )
        )
        {
            if (Pstream::parRun())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                const symmTensorField& f =
                    DEpsilonPf_.boundaryField()[patchI];

                OPstream::write
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(f.begin()),
                    f.byteSize()
                );
            }
        }
    }

    forAll(DEpsilonPf_.boundaryField(), patchI)
    {
        if
        (
            isA<processorFvsPatchSymmTensorField>
            (
                DEpsilonPf_.boundaryField()[patchI]
            )
        )
        {
            if (Pstream::parRun())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                ngbProcDEpsilonP.set
                (
                    patchI,
                    new symmTensorField
                    (
                        DEpsilonPf_.boundaryField()[patchI].size(),
                        symmTensor::zero
                    )
                );

                IPstream::read
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>
                    (
                        ngbProcDEpsilonP[patchI].begin()
                    ),
                    ngbProcDEpsilonP[patchI].byteSize()
                );
            }
        }
    }

    // DSigmaY
    forAll(DSigmaYf_.boundaryField(), patchI)
    {
        if
        (
            isA<processorFvsPatchScalarField>
            (
                DSigmaYf_.boundaryField()[patchI]
            )
        )
        {
            if (Pstream::parRun())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                const scalarField& f =
                    DSigmaYf_.boundaryField()[patchI];

                OPstream::write
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(f.begin()),
                    f.byteSize()
                );
            }
        }
    }

    forAll(DSigmaYf_.boundaryField(), patchI)
    {
        if
        (
            isA<processorFvsPatchScalarField>
            (
                DSigmaYf_.boundaryField()[patchI]
            )
        )
        {
            if (Pstream::parRun())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                ngbProcDSigmaY.set
                (
                    patchI,
                    new scalarField
                    (
                        DSigmaYf_.boundaryField()[patchI].size(),
                        0
                    )
                );

                IPstream::read
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>
                    (
                        ngbProcDSigmaY[patchI].begin()
                    ),
                    ngbProcDSigmaY[patchI].byteSize()
                );
            }
        }
    }
}
