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

#include "layerSmooth.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "attachDetach.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::layerSmooth::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston motion

    if
    (
/*
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
     || topoChanger_.size() > 0
*/
        pointZones().size() > 0
     && faceZones().size() > 0
     && topoChanger_.size() > 0
    )
    {
        Info<< "Time = " << engTime().theta() << endl;
        Info<< "void Foam::layerSmooth::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        setVirtualPositions();
        checkAndCalculate();

        Info << "Point zones found = " << pointZones().size() << endl;
        Info << "Face zones found = " << faceZones().size() << endl;
        Info << "Cell zones found = " << cellZones().size() << endl;

        return;
    }
    else
    {
        pointZones().setSize(0);
        faceZones().setSize(0);
        cellZones().setSize(0);
        topoChanger_.setSize(0);
    }

    if
    (
        engTime().engineDict().found("zOffsetGambit")
     && engTime().engineDict().found("zDisplGambit")
    )
    {
        Info << "Assembling the cylinder mesh" << endl;

        scalar zOffset
        (
            readScalar(engTime().engineDict().lookup("zOffsetGambit"))
        );

        scalar zDispl
        (
            readScalar(engTime().engineDict().lookup("zDisplGambit"))
        );

        const pointField& ap = allPoints();

        pointField pDispl = allPoints();

        forAll(ap, pointI)
        {
            const point p = ap[pointI];

            if (p.z() >= zOffset)
            {
                pDispl[pointI].z() -= zDispl;
            }
        }


        movePoints(pDispl);
        write();
        resetMotion();

        Info << "Cylinder mesh assembled" << endl;
    }



    Info << "checkAndCalculate()" << endl;
    checkAndCalculate();

    Info<< "Time = " << engTime().theta() << endl
        << "Adding zones to the engine mesh" << endl;

/*
    Point zones
    1) Piston points
*/

    DynamicList<pointZone*> pz;

/*
    Face zones
    1) Piston layer faces

*/
    DynamicList<faceZone*> fz;

/*
*/

    DynamicList<cellZone*> cz;

    label nPointZones = 0;
    label nFaceZones = 0;
    label nCellZones = 0;

/* 

    Adding the following faces zones:
    1:  pistonLayerFaces

    Adding the following point zones:
    1:  pistonPoints

*/

#   include "addPistonFacesPointZonesLayerSmooth.H"


#   include "addAttachDetachFacesLayerSmooth.H"

    Info<< "Adding " << nPointZones << " point zones, "
        << nFaceZones << " face zones and "
        << nCellZones << " cell zones" << endl;

    pz.setSize(nPointZones);
    Info << "setSize pz" << endl;
    fz.setSize(nFaceZones);
    Info << "setSize fz" << endl;
    cz.setSize(nCellZones);
    Info << "setSize cz" << endl;

    addZones(pz, fz, cz);

#   include "addMeshModifiersLayerSmooth.H"

    // Calculating the virtual positions of piston and valves

    setVirtualPositions();

    // Write mesh
    write();

    Info << "virtualPistonPosition = " << virtualPistonPosition() << endl;
    Info << "piston position = " << pistonPosition() << endl;
}


// ************************************************************************* //
