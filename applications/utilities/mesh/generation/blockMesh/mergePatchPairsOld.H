
    // Read in a list of dictionaries for the merge patch pairs
    if (meshDict.found("mergePatchPairs"))
    {
        List<Pair<word> > mergePatchPairs
        (
            meshDict.lookup("mergePatchPairs")
        );

        if (mergePatchPairs.size() > 0)
        {
            // Create and add point and face zones and mesh modifiers
            List<pointZone*> pz(mergePatchPairs.size());
            List<faceZone*> fz(3*mergePatchPairs.size());
            List<cellZone*> cz(0);

            forAll (mergePatchPairs, pairI)
            {
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                pz[pairI] = new pointZone
                (
                    mergeName + "CutPointZone",
                    labelList(0),
                    0,
                    mesh.pointZones()
                );

                // Master patch
                const word masterPatchName(mergePatchPairs[pairI].first());
                const polyPatch& masterPatch =
                    mesh.boundaryMesh()
                    [
                        mesh.boundaryMesh().findPatchID(masterPatchName)
                    ];

                labelList isf(masterPatch.size());

                forAll (isf, i)
                {
                    isf[i] = masterPatch.start() + i;
                }

                fz[3*pairI] = new faceZone
                (
                    mergeName + "MasterZone",
                    isf,
                    boolList(masterPatch.size(), false),
                    0,
                    mesh.faceZones()
                );

                // Slave patch
                const word slavePatchName(mergePatchPairs[pairI].second());
                const polyPatch& slavePatch =
                    mesh.boundaryMesh()
                    [
                        mesh.boundaryMesh().findPatchID(slavePatchName)
                    ];

                labelList osf(slavePatch.size());

                forAll (osf, i)
                {
                    osf[i] = slavePatch.start() + i;
                }

                fz[3*pairI + 1] = new faceZone
                (
                    mergeName + "SlaveZone",
                    osf,
                    boolList(slavePatch.size(), false),
                    1,
                    mesh.faceZones()
                );

                // Add empty zone for cut faces
                fz[3*pairI + 2] = new faceZone
                (
                    mergeName + "CutFaceZone",
                    labelList(0),
                    boolList(0, false),
                    2,
                    mesh.faceZones()
                );
            }  // end of all merge pairs

            Info << "Adding point and face zones" << endl;
            mesh.addZones(pz, fz, cz);


            Info << "Creating topo change" << endl;
            polyTopoChanger attacher(mesh);
            attacher.setSize(mergePatchPairs.size());

            forAll (mergePatchPairs, pairI)
            {
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                // Add the sliding interface mesh modifier
                attacher.set
                (
                    pairI,
                    new slidingInterface
                    (
                        "couple" + name(pairI),
                        pairI,
                        attacher,
                        mergeName + "MasterZone",
                        mergeName + "SlaveZone",
                        mergeName + "CutPointZone",
                        mergeName + "CutFaceZone",
                        mergePatchPairs[pairI].first(),
                        mergePatchPairs[pairI].second(),
                        slidingInterface::INTEGRAL, // always integral
                        intersection::VISIBLE
                    )
                );
            }

            attacher.changeMesh();

            // Clean the mesh after attach
            labelList patchSizes(mesh.boundaryMesh().size());
            labelList patchStarts(mesh.boundaryMesh().size());

            forAll (mesh.boundaryMesh(), patchI)
            {
                patchSizes[patchI] = mesh.boundaryMesh()[patchI].size();
                patchStarts[patchI] = mesh.boundaryMesh()[patchI].start();
            }

            mesh.resetPrimitives
            (
                xferCopy<pointField>(mesh.points()),
                xferCopy<faceList>(mesh.faces()),
                xferCopy<labelList>(mesh.faceOwner()),
                xferCopy<labelList>(mesh.faceNeighbour()),
                patchSizes,
                patchStarts
            );

            mesh.setInstance(runTime.constant());
            mesh.removeZones();
        }
    }
    else
    {
        Info<< nl << "There are no merge patch pairs" << endl;
    }


    // Set any cellZones (note: cell labelling unaffected by above
    // mergePatchPairs)

    label nZones = blocks.numZonedBlocks();

    if (nZones > 0)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<dynamicLabelList > zoneCells(nZones);

        // Running cell counter
        label cellI = 0;

        // Largest zone so far
        label freeZoneI = 0;

        forAll(blocks, blockI)
        {
            const block& b = blocks[blockI];
            const labelListList& blockCells = b.cells();
            const word& zoneName = b.blockDef().zoneName();

            if (zoneName.size())
            {
                HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

                label zoneI;

                if (iter == zoneMap.end())
                {
                    zoneI = freeZoneI++;

                    Info<< "    " << zoneI << '\t' << zoneName << endl;

                    zoneMap.insert(zoneName, zoneI);
                }
                else
                {
                    zoneI = iter();
                }

                forAll(blockCells, i)
                {
                    zoneCells[zoneI].append(cellI++);
                }
            }
            else
            {
                cellI += b.cells().size();
            }
        }


        List<cellZone*> cz(zoneMap.size());

        Info<< nl << "Writing cell zones as cellSets" << endl;

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI]= new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                zoneI,
                mesh.cellZones()
            );

            // Write as cellSet for ease of processing
            cellSet cset
            (
                mesh,
                iter.key(),
                labelHashSet(zoneCells[zoneI].shrink())
            );
            cset.write();
        }

        mesh.pointZones().setSize(0);
        mesh.faceZones().setSize(0);
        mesh.cellZones().setSize(0);
        mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << nl << "Writing polyMesh" << endl;
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
