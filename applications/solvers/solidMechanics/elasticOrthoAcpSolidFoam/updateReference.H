{
    // Check which region needs reference
    needRef = boolList(mesh.regions().nRegions(), true);

    forAll(U.boundaryField(), patchI)
    {
        // philipc - this used to set a reference on
        // processors which did not have a patch that fixesValue
        // so processor in the centre of the domain all had
        // a referece set causing stress peaks and resulting
        // in an incorrect solution
        // so a quick fix is to not set a reference on regions
        // with a processor boundary
        //if (U.boundaryField()[patchI].fixesValue())
        if
        (
            U.boundaryField()[patchI].fixesValue()
            ||
            mesh.boundaryMesh()[patchI].type()
            == processorPolyPatch::typeName
        )
        {
            const unallocLabelList& curFaceCells =
                mesh.boundary()[patchI].faceCells();

            forAll(curFaceCells, fI)
            {
                needRef[mesh.regions()[curFaceCells[fI]]] = false;
            }
        }
    }

    // Find reference cell for each region
    refCells = labelList(mesh.regions().nRegions(), -1);

    forAll(refCells, regI)
    {
        refCells[regI] = findIndex(mesh.regions(), regI);

        if (refCells[regI] == -1)
        {
            FatalErrorIn(args.executable())
                << "Can't find reference cell for region: "
                << regI << abort(FatalError);
        }
    }

    // Find dead cells
    labelHashSet deadCellsSet;
    const labelList& regs = mesh.regions();

    forAll(regs, cellI)
    {
        if (mesh.nCellsInRegion(regs[cellI]) < 0.01*mesh.nCells())
        {
            deadCellsSet.insert(cellI);
        }
    }

    deadCells = deadCellsSet.toc();
}
