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

#include "simpleTwoStroke.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleTwoStroke::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston motion

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "Time = " << engTime().theta() << endl;
        Info<< "void simpleTwoStroke::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void simpleTwoStroke::addZonesAndModifiers()"
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

    List<pointZone*> pz(3);
    List<faceZone*> fz(4);
    List<cellZone*> cz(1);

    label nPointZones = 0;
    label nFaceZones = 0;
    label nCellZones = 0;

    // Add the piston zone
    if (piston().patchID().active())
    {

        // Piston position

        Info << "Adding face zone for piston layer addition/removal" << endl;

        label pistonPatchID = piston().patchID().index();

        scalar zPist =
            max(boundary()[pistonPatchID].patch().localPoints()).z();

        scalar zPistV = zPist + offSet();

        labelList zone1(faceCentres().size());
        boolList flipZone1(faceCentres().size(), false);
        label nZoneFaces1 = 0;

        bool foundAtLeastOne = false;
        scalar zHigher = GREAT;
        scalar zLower = GREAT;
        scalar dh = GREAT;
        scalar dl = GREAT;

        forAll (faceCentres(), faceI)
        {
            // The points have to be in the cylinder and not in the ports....

            scalar zc = faceCentres()[faceI].z();

            scalar xc = faceCentres()[faceI].x();
            scalar yc = faceCentres()[faceI].y();

            vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);
            scalar dd = n & vector(0,0,1);

            if(sqrt(sqr(xc)+sqr(yc)) <  0.5 * engTime().bore().value())
            {
                if (dd > 0.1)
                {
                    if (zPistV - zc > 0 && zPistV - zc < dl)
                    {
                        zLower = zc;
                        dl = zPistV - zc;
                    }

                    if (zc - zPistV > 0 && zc - zPistV < dh)
                    {
                        zHigher = zc;
                        dh = zc - zHigher;
                    }

                    if
                    (
                        zc > zPistV - delta()
                        && zc < zPistV + delta()
                    )
                    {
                        foundAtLeastOne = true;
                        if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                        {
                            flipZone1[nZoneFaces1] = true;
                        }

                        zone1[nZoneFaces1] = faceI;
                        nZoneFaces1++;
                    }
                }
            }
        }

        // if no cut was found use the layer above
        if (!foundAtLeastOne)
        {
            zPistV = zHigher;

            forAll (faceCentres(), faceI)
            {
                scalar zc = faceCentres()[faceI].z();

                scalar xc = faceCentres()[faceI].x();
                scalar yc = faceCentres()[faceI].y();

                vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);
                scalar dd = n & vector(0,0,1);

                if(sqrt(sqr(xc)+sqr(yc)) <  0.5 * engTime().bore().value())
                {
                    if (dd > 0.1)
                    {

                        if
                        (
                            zc > zPistV - delta()
                            && zc < zPistV + delta()
                        )
                        {
                            if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                            {
                                flipZone1[nZoneFaces1] = true;
                            }

                            zone1[nZoneFaces1] = faceI;
                            nZoneFaces1++;
                        }
                    }
                }

            }
        }

        zone1.setSize(nZoneFaces1);
        flipZone1.setSize(nZoneFaces1);

        fz[nFaceZones] =
            new faceZone
            (
                "pistonLayerFaces",
                zone1,
                flipZone1,
                nFaceZones,
                faceZones()
            );

        nFaceZones++;

        // Points which don't move (= cylinder head)
        DynamicList<label> headPoints(nPoints() / 10);

        label nHeadPoints = 0;
        forAll (points(), pointI)
        {
            scalar zCoord = points()[pointI].z();

            if (zCoord > deckHeight() - delta())
            {
                headPoints.append(pointI);
                nHeadPoints++;
            }
        }

        Info << "Number of head points = " << nHeadPoints << endl;
        pz[nPointZones] =
            new pointZone
            (
                "headPoints",
                headPoints.shrink(),
                nPointZones,
                pointZones()
            );

        nPointZones++;

    }

    //  Sliding interface for scavenging ports

    if (foundScavPorts())
    {
        // Inner slider

        const polyPatch& innerScav =
            boundaryMesh()[boundaryMesh().findPatchID(scavInCylPatchName_)];

        labelList isf(innerScav.size());

        forAll (isf, i)
        {
            isf[i] = innerScav.start() + i;
        }

        fz[nFaceZones] = new faceZone
        (
            scavInCylPatchName_ + "Zone",
            isf,
            boolList(innerScav.size(), false),
            nFaceZones,
            faceZones()
        );

        nFaceZones++;

        // Outer slider

        const polyPatch& outerScav =
            boundaryMesh()[boundaryMesh().findPatchID(scavInPortPatchName_)];

        labelList osf(outerScav.size());

        forAll (osf, i)
        {
            osf[i] = outerScav.start() + i;
        }

        fz[nFaceZones] = new faceZone
        (
            scavInPortPatchName_ + "Zone",
            osf,
            boolList(outerScav.size(), false),
            nFaceZones,
            faceZones()
        );

        nFaceZones++;

        fz[nFaceZones] = new faceZone
        (
            "cutFaceZone",
            labelList(0),
            boolList(0, false),
            nFaceZones,
            faceZones()
        );

        nFaceZones++;

        Info << "cut p" << endl;

        pz[nPointZones] = new pointZone
        (
            "cutPointZone",
            labelList(0),
            nPointZones,
            pointZones()
        );

        nPointZones++;
    }


    {
        labelList movingCells(nCells());
        label nMovingCells = 0;

        scalar pistonX = piston().cs().origin().x();
        scalar pistonY = piston().cs().origin().y();

        forAll(cellCentres(),cellI)
        {
            const vector& v = cellCentres()[cellI];

            if
            (
                sqrt(sqr(v.x()-pistonX)+sqr(v.y()-pistonY))
             < 0.5*engTime().bore().value()
            )
            {
                movingCells[nMovingCells] = cellI;
                nMovingCells++;
            }
        }

        movingCells.setSize(nMovingCells);
        Info<< "Number of cells in the moving region: "
            << nMovingCells << endl;

        cz[nCellZones] = new cellZone
        (
            "movingCells",
            movingCells,
            nCellZones,
            cellZones()
        );

        nCellZones++;
    }


    Info<< "Adding " << nPointZones << " point, "
        << nFaceZones << " face zones and "
        << nCellZones << " cell zones" << endl;

    pz.setSize(nPointZones);
    fz.setSize(nFaceZones);
    cz.setSize(nCellZones);
    addZones(pz, fz, cz);

    List<polyMeshModifier*> tm(2);
    label nMods = 0;

    // Add piston layer addition
    if (piston().patchID().active())
    {

        topoChanger_.setSize(topoChanger_.size() + 1);

        topoChanger_.set
        (
            nMods,
            new layerAdditionRemoval
            (
                "pistonLayer",
                nMods,
                topoChanger_,
                "pistonLayerFaces",
                piston().minLayer(),
                piston().maxLayer()
            )
        );
        nMods++;
    }


    if(foundScavPorts())
    {
        topoChanger_.setSize(topoChanger_.size() + 1);

        topoChanger_.set
        (
            nMods,
            new slidingInterface
            (
                "valveSlider",
                nMods,
                topoChanger_,
                scavInCylPatchName_ + "Zone",
                scavInPortPatchName_ + "Zone",
                "cutPointZone",
                "cutFaceZone",
                scavInCylPatchName_,
                scavInPortPatchName_,
                slidingInterface::INTEGRAL,
//              slidingInterface::PARTIAL,
                true,                          // Attach-detach action
                intersection::VISIBLE         // Projection algorithm
            )
        );
        nMods++;
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
