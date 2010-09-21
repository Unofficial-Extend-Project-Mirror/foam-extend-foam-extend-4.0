/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    blockMesh

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
    @a constant/polyMesh/blockMeshDict
    (or @a constant/\<region\>/polyMesh/blockMeshDict).

Usage

    - blockMesh [OPTION]

    @param -blockTopology \n
    Write the topology as a set of edges in OBJ format.

    @param -region \<name\> \n
    Specify an alternative mesh region.

    @param -dict \<dictionary\> \n
    Specify an alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "mapPolyMesh.H"
#include "polyTopoChanger.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("blockTopology", "");
    argList::validOptions.insert("dict", "dictionary");
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    const word dictName("blockMeshDict");

    word regionName;
    fileName polyMeshDir;

    if (args.optionFound("region"))
    {
        // constant/<region>/polyMesh/blockMeshDict
        regionName  = args.option("region");
        polyMeshDir = regionName/polyMesh::meshSubDir;

        Info<< nl << "Generating mesh for region " << regionName << endl;
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        regionName  = polyMesh::defaultRegion;
        polyMeshDir = polyMesh::meshSubDir;
    }

    autoPtr<IOobject> meshDictIoPtr;

    if (args.optionFound("dict"))
    {
        fileName dictPath(args.option("dict"));

        meshDictIoPtr.set
        (
            new IOobject
            (
                (
                    isDir(dictPath)
                  ? dictPath/dictName
                  : dictPath
                ),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
    }
    else
    {
        meshDictIoPtr.set
        (
            new IOobject
            (
                dictName,
                runTime.constant(),
                polyMeshDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
    }

    if (!meshDictIoPtr->headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot open mesh description file\n    "
            << meshDictIoPtr->objectPath()
            << nl
            << exit(FatalError);
    }

    Info<< nl << "Creating block mesh from\n    "
        << meshDictIoPtr->objectPath() << nl << endl;

    IOdictionary meshDict(meshDictIoPtr());
    blockMesh blocks(meshDict);


    if (args.optionFound("blockTopology"))
    {
        // Write mesh as edges.
        {
            fileName objMeshFile("blockTopology.obj");

            OFstream str(runTime.path()/objMeshFile);

            Info<< nl << "Dumping block structure as Lightwave obj format"
                << " to " << objMeshFile << endl;

            blocks.writeTopology(str);
        }

        // Write centres of blocks
        {
            fileName objCcFile("blockCentres.obj");

            OFstream str(runTime.path()/objCcFile);

            Info<< nl << "Dumping block centres as Lightwave obj format"
                << " to " << objCcFile << endl;

            const polyMesh& topo = blocks.topology();

            const pointField& cellCentres = topo.cellCentres();

            forAll(cellCentres, cellI)
            {
                //point cc = b.blockShape().centre(b.points());
                const point& cc = cellCentres[cellI];

                str << "v " << cc.x() << ' ' << cc.y() << ' ' << cc.z() << nl;
            }
        }

        Info<< nl << "end" << endl;

        return 0;
    }



    Info<< nl << "Creating mesh from block mesh" << endl;

    wordList patchNames = blocks.patchNames();
    wordList patchTypes = blocks.patchTypes();
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    wordList patchPhysicalTypes = blocks.patchPhysicalTypes();

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMeshDir,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        xferCopy(blocks.points()),
        blocks.cells(),
        blocks.patches(),
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );


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
                        slidingInterface::INTEGRAL,     // always integral
                        false,                          // attach-detach action
                        intersection::VISIBLE
                    )
                );
            }

            attacher.changeMesh();
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
        List<DynamicList<label> > zoneCells(nZones);

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
