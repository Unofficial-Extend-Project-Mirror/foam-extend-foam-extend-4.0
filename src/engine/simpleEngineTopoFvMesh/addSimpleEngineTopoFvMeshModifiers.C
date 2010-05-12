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

#include "simpleEngineTopoFvMesh.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "attachDetach.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleEngineTopoFvMesh::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston and valve motion

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void Foam::simpleEngineTopoFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void simpleEngineTopoFvMesh::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        return;
    }

    Info<< "Time = " << engineTime_.theta() << endl
        << "Adding zones to the engine mesh" << endl;

    List<pointZone*> pz(nValves());
    List<faceZone*> fz(6*nValves() + 1);
    List<cellZone*> cz(0);

    label nPointZones = 0;
    label nFaceZones = 0;

    for (label valveI = 0; valveI < nValves(); valveI++)
    {
        // If both sides of the interface exist, add sliding interface
        // for a valve
        if
        (
            valves_[valveI].curtainInCylinderPatchID().active()
         && valves_[valveI].curtainInPortPatchID().active()
        )
        {
            Info<< "Adding sliding interface zones for curtain of valve "
                << valveI + 1 << endl;

            pz[nPointZones] =
                new pointZone
                (
                    "cutPointsV" + Foam::name(valveI + 1),
                    labelList(0),
                    nPointZones,
                    pointZones()
                );
            nPointZones++;

            const polyPatch& cylCurtain =
                boundaryMesh()
                    [valves_[valveI].curtainInCylinderPatchID().index()];

            labelList cylCurtainLabels(cylCurtain.size(), cylCurtain.start());

            forAll (cylCurtainLabels, i)
            {
                cylCurtainLabels[i] += i;
            }

            fz[nFaceZones] =
                new faceZone
                (
                    "curtainCylZoneV" + Foam::name(valveI + 1),
                    cylCurtainLabels,
                    boolList(cylCurtainLabels.size(), false),
                    nFaceZones,
                    faceZones()
                );
            nFaceZones++;

            const polyPatch& portCurtain =
                boundaryMesh()
                    [valves_[valveI].curtainInPortPatchID().index()];

            labelList portCurtainLabels
            (
                portCurtain.size(),
                portCurtain.start()
            );

            forAll (portCurtainLabels, i)
            {
                portCurtainLabels[i] += i;
            }

            fz[nFaceZones] =
                new faceZone
                (
                    "curtainPortZoneV" + Foam::name(valveI + 1),
                    portCurtainLabels,
                    boolList(portCurtainLabels.size(), false),
                    nFaceZones,
                    faceZones()
                );
            nFaceZones++;

            // Add empty zone for cut faces
            fz[nFaceZones] =
                new faceZone
                (
                    "cutFaceZoneV" + Foam::name(valveI + 1),
                    labelList(0),
                    boolList(0, false),
                    nFaceZones,
                    faceZones()
                );
            nFaceZones++;

            // Create a detach zone
            if
            (
                valves_[valveI].detachInCylinderPatchID().active()
             && valves_[valveI].detachInPortPatchID().active()
             && valves_[valveI].detachFaces().size() > 0
            )
            {
                Info<< "Adding detach boundary for valve "
                    << valveI + 1 << endl;

                const vectorField& areas = Sf().internalField();

                const labelList& df = valves_[valveI].detachFaces();

                boolList flip(df.size(), false);

                const vector& pistonAxis = piston().cs().axis();

                forAll (df, dfI)
                {
                    if (isInternalFace(df[dfI]))
                    {
                        if ((areas[df[dfI]] & pistonAxis) > 0)
                        {
                            flip[dfI] = true;
                        }
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "void simpleEngineTopoFvMesh::"
                            "addZonesAndModifiers()"
                        )   << "found boundary face in valve detach definition"
                            << " for valve " << valveI + 1
                            << ".  This is not allowed.  Detach faces: "
                            << df << " nInternalFaces: " << nInternalFaces()
                            << abort(FatalError);
                    }
                }

                // Add detach face zone
                fz[nFaceZones] =
                    new faceZone
                    (
                        "detachFaceZoneV" + Foam::name(valveI + 1),
                        df,
                        flip,
                        nFaceZones,
                        faceZones()
                    );
                nFaceZones++;
            }
        }
        else
        {
            Info << "No valve curtain for valve " << valveI + 1 << endl;
        }

        // Make a zone for layer addition at the top of the valve
        if (valves_[valveI].poppetPatchID().active())
        {
            Info<< "Adding poppet layer addition zone for valve "
                << valveI + 1 << endl;

            const polyPatch& poppetPatch =
                boundaryMesh()
                    [valves_[valveI].poppetPatchID().index()];

            labelList poppetPatchLabels
            (
                poppetPatch.size(),
                poppetPatch.start()
            );

            forAll (poppetPatchLabels, i)
            {
                poppetPatchLabels[i] += i;
            }

            fz[nFaceZones] =
                new faceZone
                (
                    "poppetZoneV" + Foam::name(valveI + 1),
                    poppetPatchLabels,
                    boolList(poppetPatchLabels.size(), true),
                    nFaceZones,
                    faceZones()
                );
            nFaceZones++;
        }
        else
        {
            Info << "No poppet layer addition zone for valve "
                << valveI + 1 << endl;
        }

        if (valves_[valveI].bottomPatchID().active())
        {
            Info<< "Adding bottom layer addition zone for valve "
                << valveI + 1 << endl;

            const polyPatch& bottomPatch =
                boundaryMesh()
                    [valves_[valveI].bottomPatchID().index()];

            labelList bottomPatchLabels
            (
                bottomPatch.size(),
                bottomPatch.start()
            );

            forAll (bottomPatchLabels, i)
            {
                bottomPatchLabels[i] += i;
            }

            fz[nFaceZones] =
                new faceZone
                (
                    "bottomZoneV" + Foam::name(valveI + 1),
                    bottomPatchLabels,
                    boolList(bottomPatchLabels.size(), true),
                    nFaceZones,
                    faceZones()
                );
            nFaceZones++;
        }
        else
        {
            Info << "No bottom layer addition zone for valve "
                << valveI + 1 << endl;
        }
    }

    // Add the piston zone
    if (piston().patchID().active())
    {
        Info << "Adding layer addition zone for piston" << endl;

        const polyPatch& pistonPatch =
            boundaryMesh()[piston().patchID().index()];

        labelList pistonPatchLabels(pistonPatch.size(), pistonPatch.start());

        forAll (pistonPatchLabels, i)
        {
            pistonPatchLabels[i] += i;
        }

        fz[nFaceZones] =
            new faceZone
            (
                "pistonZone",
                pistonPatchLabels,
                boolList(pistonPatchLabels.size(), true),
                nFaceZones,
                faceZones()
            );
        nFaceZones++;
    }

    Info<< "Adding " << nPointZones << " point and "
        << nFaceZones << " face zones" << endl;

    pz.setSize(nPointZones);
    fz.setSize(nFaceZones);
    addZones(pz, fz, cz);

    topoChanger_.setSize(4*nValves() + 1);
    label nMods = 0;

    for (label valveI = 0; valveI < nValves(); valveI++)
    {
        // Add valve curtain sliding interface
        if
        (
            valves_[valveI].curtainInCylinderPatchID().active()
         && valves_[valveI].curtainInPortPatchID().active()
        )
        {
            topoChanger_.set
            (
                nMods,
                new slidingInterface
                (
                    "valveSlider" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "curtainPortZoneV" + Foam::name(valveI + 1),
                    "curtainCylZoneV" + Foam::name(valveI + 1),
                    "cutPointsV" + Foam::name(valveI + 1),
                    "cutFaceZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].curtainInPortPatchID().name(),
                    valves_[valveI].curtainInCylinderPatchID().name(),
                    slidingInterface::INTEGRAL, // always integral
                    true,  // attach-detach action
                    intersection::VISIBLE
                )
            );
            nMods++;
        }

        // Add attach-detach for valve
        if
        (
            valves_[valveI].detachInCylinderPatchID().active()
         && valves_[valveI].detachInPortPatchID().active()
         && valves_[valveI].detachFaces().size() > 0
        )
        {
            topoChanger_.set
            (
                nMods,
                new attachDetach
                (
                    "valveAttachDetach" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "detachFaceZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].detachInCylinderPatchID().name(),
                    valves_[valveI].detachInPortPatchID().name(),
                    scalarField(0),
                    true                // Manual triggering
                )
            );
            nMods++;
        }
        // Add valve poppet layer addition
        if (valves_[valveI].poppetPatchID().active())
        {
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    "valvePoppetLayer" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "poppetZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].minTopLayer(),
                    valves_[valveI].maxTopLayer()
                )
            );
            nMods++;
        }

        // Add valve bottom layer addition
        if (valves_[valveI].bottomPatchID().active())
        {
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    "valveBottomLayer" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "bottomZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].minBottomLayer(),
                    valves_[valveI].maxBottomLayer()
                )
            );
            nMods++;
        }
    }

    // Add piston layer addition
    if (piston().patchID().active())
    {
        topoChanger_.set
        (
            nMods,
            new layerAdditionRemoval
            (
                "pistonLayer",
                nMods,
                topoChanger_,
                "pistonZone",
                piston().minLayer(),
                piston().maxLayer()
            )
        );
        nMods++;
    }

    Info << "Adding " << nMods << " topology modifiers" << endl;
    topoChanger_.setSize(nMods);

    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


// ************************************************************************* //
