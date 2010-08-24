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

#include "layerAR.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::layerAR::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston motion

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void layerAR::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void layerAR::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        setVirtualPistonPosition();
        checkAndCalculate();

        return;
    }

    checkAndCalculate();
    
    Info<< "Time = " << engTime().theta() << endl
        << "Adding zones to the engine mesh" << endl;

    //fz = 1: faces where layer are added/removed
    //pz = 2: points below the virtual piston faces and head points

    List<pointZone*> pz(2);
    List<faceZone*> fz(1);
    List<cellZone*> cz(0);

    label nPointZones = 0;
    label nFaceZones = 0;

    // Add the piston zone
    if (piston().patchID().active() && offSet() > SMALL)
    {

        // Piston position
        
        label pistonPatchID = piston().patchID().index();
        
        scalar zPist = max(boundary()[pistonPatchID].patch().localPoints()).z();
        
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
            scalar zc = faceCentres()[faceI].z();
            vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);
            scalar dd = n & vector(0,0,1);

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

        // if no cut was found use the layer above
        if (!foundAtLeastOne)
        {
                        
            zPistV = zHigher;

            forAll (faceCentres(), faceI)
            {
                scalar zc = faceCentres()[faceI].z();
                vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);
                scalar dd = n & vector(0,0,1);
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

        zone1.setSize(nZoneFaces1);
        flipZone1.setSize(nZoneFaces1);
    
        fz[nFaceZones]=
            new faceZone
            (
                "pistonLayerFaces",
                zone1,
                flipZone1,
                nFaceZones,
                faceZones()
            );
        
        nFaceZones++;


        // Construct point zones

            
        // Points which don't move (= cylinder head)
        DynamicList<label> headPoints(nPoints() / 10);

        // Points below the piston which moves with the piston displacement
        DynamicList<label> pistonPoints(nPoints() / 10);
        
        label nHeadPoints = 0;
            
        forAll (points(), pointI)
        {
            scalar zCoord = points()[pointI].z();

            if (zCoord > deckHeight() - delta())
            {
                headPoints.append(pointI);
                nHeadPoints++; 
            }
            else if (zCoord < zPistV + delta())
            {
                pistonPoints.append(pointI);
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

        pz[nPointZones] =
            new pointZone
            (
                "pistonPoints",
                pistonPoints.shrink(),
                nPointZones,
                pointZones()
            );

        nPointZones++;

    }
    else if(piston().patchID().active() && offSet() <= SMALL)
    {
        label pistonPatchID = piston().patchID().index();
        
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
                "pistonLayerFaces",
                pistonPatchLabels,
                boolList(pistonPatchLabels.size(), true),
                nFaceZones,
                faceZones()
            );
        nFaceZones++;
        // Construct point zones

        scalar zPistV = max(boundary()[pistonPatchID].patch().localPoints()).z();
            
        // Points which don't move (= cylinder head)
        DynamicList<label> headPoints(nPoints() / 10);

        // Points below the piston which moves with the piston displacement
        DynamicList<label> pistonPoints(nPoints() / 10);
        
        label nHeadPoints = 0;
            
        forAll (points(), pointI)
        {
            scalar zCoord = points()[pointI].z();

            if (zCoord > deckHeight() - delta())
            {
                headPoints.append(pointI);
                nHeadPoints++; 
            }
            else if (zCoord < zPistV + delta())
            {
                pistonPoints.append(pointI);
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

        pz[nPointZones] =
            new pointZone
            (
                "pistonPoints",
                pistonPoints.shrink(),
                nPointZones,
                pointZones()
            );

        nPointZones++;

    }


    Info<< "Adding " << nPointZones << " point and "
        << nFaceZones << " face zones" << endl;

    pz.setSize(nPointZones);
    fz.setSize(nFaceZones);
    addZones(pz, fz, cz);

    List<polyMeshModifier*> tm(1);
    label nMods = 0;

    // Add piston layer addition
    Info << "Adding Layer Addition/Removal Mesh Modifier" << endl;

    if (piston().patchID().active())
    {
        topoChanger_.setSize(1);
        topoChanger_.set
        (
            0,
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
    }

    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();

    // Calculating the virtual piston position
    setVirtualPistonPosition();

    Info << "virtualPistonPosition = " << virtualPistonPosition() << endl;
    Info << "piston position = " << pistonPosition() << endl;
}


// ************************************************************************* //
