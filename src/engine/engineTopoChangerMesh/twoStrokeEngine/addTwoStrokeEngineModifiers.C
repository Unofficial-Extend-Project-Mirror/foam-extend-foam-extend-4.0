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

\*---------------------------------------------------------------------------*/

#include "twoStrokeEngine.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "SortableList.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoStrokeEngine::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston motion
    if (faceZones().size() > 0)
    {
        Info<< "Time = " << engTime().theta() << endl;
        Info<< "void twoStrokeEngine::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void twoStrokeEngine::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        setVirtualPistonPosition();
        checkAndCalculate();

        return;
    }

    Info << "checkAndCalculate()" << endl;
    checkAndCalculate();

    Info<< "Time = " << engTime().theta() << endl
        << "Adding zones to the engine mesh" << endl;


    //fz = 4: virtual piston, outSidePort, insidePort, cutFaceZone
    //pz = 2: piston points, cutPointZone
    //cz = 1: moving mask

    label nPorts = scavInCylPatches_.size();


    DynamicList<pointZone*> pz(2 + nPorts);
    DynamicList<faceZone*> fz(3*nPorts + 1);

    // Added piston cells and head cells
    DynamicList<cellZone*> cz(3);

    label nPointZones = 0;
    label nFaceZones = 0;
    label nCellZones = 0;

    Info << "Adding piston layer faces" << endl;

#   include "addPistonLayerHrv.H"

//  adding head points (does not move)

    {
//         pointSet headPointSet(*this, headPointsSetName_);

//         pz[nPointZones] =
//             new pointZone
//             (
//                 "headPoints",
//                 headPointSet.toc(),
//                 nPointZones,
//                 pointZones()
//             );

//         nPointZones++;

        cellSet headCellSet(*this, headCellsSetName_);

        cz[nCellZones] =
            new cellZone
            (
                "headCells",
                headCellSet.toc(),
                nCellZones,
                cellZones()
            );

        nCellZones++;
    }

//  Sliding interface for scavenging ports

    if(nPorts > 0)
    {
        forAll(scavInCylPatches_, patchi)
        {
            // Inner slider
            const polyPatch& innerScav =
                boundaryMesh()
                [
                    boundaryMesh().findPatchID(scavInCylPatches_[patchi])
                ];

            labelList isf(innerScav.size());

            forAll (isf, i)
            {
                isf[i] = innerScav.start() + i;
            }

            fz[nFaceZones] = new faceZone
            (
                scavInCylPatches_[patchi] + "Zone" + Foam::name(patchi + 1),
                isf,
                boolList(innerScav.size(), false),
                nFaceZones,
                faceZones()
            );

            nFaceZones++;

            // Outer slider

            const polyPatch& outerScav =
                boundaryMesh()
                [
                    boundaryMesh().findPatchID(scavInPortPatches_[patchi])
                ];

            labelList osf(outerScav.size());

            forAll (osf, i)
            {
                osf[i] = outerScav.start() + i;
            }

            fz[nFaceZones] = new faceZone
            (
                scavInPortPatches_[patchi] + "Zone" + Foam::name(patchi + 1),
                osf,
                boolList(outerScav.size(), false),
                nFaceZones,
                faceZones()
            );

            nFaceZones++;

            fz[nFaceZones] = new faceZone
            (
                "cutFaceZone" + Foam::name(patchi + 1),
                labelList(0),
                boolList(0, false),
                nFaceZones,
                faceZones()
            );

            nFaceZones++;

            Info << "cut p" << endl;

            pz[nPointZones] = new pointZone
            (
                "cutPointZone" + Foam::name(patchi + 1),
                labelList(0),
                nPointZones,
                pointZones()
            );

            nPointZones++;
        }
    }

    Info << "moving cells" << endl;

    {
        cellSet movingCells(*this, movingCellSetName_);

        cz[nCellZones] = new cellZone
        (
            "movingCells",
            movingCells.toc(),
            nCellZones,
            cellZones()
        );

        nCellZones++;
    }

    Info << "moving cells added" << endl;

    Info<< "Adding " << nPointZones << " point, "
        << nFaceZones << " face zones and " << nCellZones
        << " cell zones" << endl;

    pz.setSize(nPointZones);
    fz.setSize(nFaceZones);
    cz.setSize(nCellZones);
    addZones(pz, fz, cz);

    List<polyMeshModifier*> tm(1 + nPorts);
    label nMods = 0;

    // Add piston layer addition

#   include "addPistonLayerAdditionRemovalMeshModifier.H"

    if(nPorts > 0)
    {
        forAll(scavInPortPatches_, i)
        {
            topoChanger_.setSize(topoChanger_.size() + 1);

            topoChanger_.set
            (
                nMods,
                new slidingInterface
                (
                    "portCylinderInterface" + Foam::name(i + 1),
                    nMods,
                    topoChanger_,
                    scavInPortPatches_[i] + "Zone" + Foam::name(i + 1),
                    scavInCylPatches_[i] + "Zone" + Foam::name(i + 1),
                    "cutPointZone" + Foam::name(i + 1),
                    "cutFaceZone" + Foam::name(i + 1),
                    scavInPortPatches_[i],
                    scavInCylPatches_[i],
                    slidingInterface::INTEGRAL,
                    true,                          // Attach-detach action
                    intersection::VISIBLE         // Projection algorithm
                )
            );

            nMods++;
        }
    }

    Info << "Adding " << nMods << " topology modifiers" << endl;

    // Calculating the virtual piston position

    setVirtualPistonPosition();


    topoChanger_.setSize(nMods);

    // Write mesh modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();

    Info << "virtualPistonPosition = " << virtualPistonPosition() << endl;
    Info << "piston position = " << pistonPosition() << endl;
}


// ************************************************************************* //
