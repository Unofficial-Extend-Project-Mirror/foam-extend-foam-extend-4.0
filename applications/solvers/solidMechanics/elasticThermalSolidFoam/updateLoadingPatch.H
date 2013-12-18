/* ------------------------------------------------------------------------- *\
||                                                                           ||
||                             Update loading patch                          ||
||                                                                           ||
\* ------------------------------------------------------------------------- */
{
    word patchName = "right";
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if(patchID == -1)
    {
        FatalErrorIn("updateLoadingPatch.H") << "Patch " << patchName << " not found." << endl
            << exit(FatalError);
    }

    //- update patch
    if
    (
        U.boundaryField()[patchID].type()
        == fixedValueFvPatchVectorField::typeName
    )
    {
        //- 1 mm every 100 seconds
        //- experimental ultimate failure occurred at 0.3 mm
        vector disp(0.00001*runTime.value(), 0, 0);
        U.boundaryField()[patchID] == disp;

        Info << "Patch " << mesh.boundary()[patchID].name() << " has been updated"
            << " to " << disp
            << endl;
    }
    else
    {
        SeriousError << "Loading Patch " << patchName << " is type "
            << U.boundaryField()[patchID].type()
            << " and is not updated!" << endl
            << exit(FatalError);
    }
}
// ************************************************************************* //



