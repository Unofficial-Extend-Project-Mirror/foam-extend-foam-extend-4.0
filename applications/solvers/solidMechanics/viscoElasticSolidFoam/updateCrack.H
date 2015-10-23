{
    forAll (DU.boundaryField(), patchI)
    {
        if
        (
            DU.boundaryField()[patchI].type()
         == cohesiveZoneIncrementalModeIIFvPatchVectorField::typeName
        )
        {
            cohesiveZoneIncrementalModeIIFvPatchVectorField& cohesivePatchDU =
                refCast<cohesiveZoneIncrementalModeIIFvPatchVectorField>
                (
                    DU.boundaryField()[patchI]
                );

            nCrackedFaces = cohesivePatchDU.updateCrack();
            Info << "nCrackedFaces is " << nCrackedFaces << endl;
        }
    }
}
