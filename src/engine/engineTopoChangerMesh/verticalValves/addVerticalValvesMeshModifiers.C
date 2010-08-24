/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2010 Tommaso Lucchini
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

Class
    verticalValves

\*---------------------------------------------------------------------------*/

#include "verticalValves.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "attachDetach.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::verticalValves::addZonesAndModifiers()
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
        Info<< "void verticalValves::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void verticalValves::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        setVirtualPositions();
        checkAndCalculate();

        Info << "Point zones found = " << pointZones().size() << endl;
        Info << "Face zones found = " << faceZones().size() << endl;
        Info << "Cell zones found = " << cellZones().size() << endl;

        return;
    }

    Info << "checkAndCalculate()" << endl;
    checkAndCalculate();

    Info<< "Time = " << engTime().theta() << endl
        << "Adding zones to the engine mesh" << endl;

/*
    Point zones
    1) Piston points
    1) Cut point zone for liner in head

    nValves*
    1) cutPointsV
    2) valveTopPoints
    3) valveBottomPoints
*/

    DynamicList<pointZone*> pz;

/*
    Face zones
    1) Piston layer faces

    nValves*
    1) valveTopLayerFaces
    2) valveBottomLayerFaces
    3) valveCurtainPort
    4) valveCurtainCyl
    5) cutFaceV
*/
    DynamicList<faceZone*> fz;

/*
    cell zones
    1) moving cells inside piston

    nValves*
    1) moving cells in the top of the valve
    2) moving cells in the bottom of the valve
*/

    DynamicList<cellZone*> cz;

    label nPointZones = 0;
    label nFaceZones = 0;
    label nCellZones = 0;

/*
    Adding the following faces zones:
    1:  pistonLayerFaces
    nV: pistonLayerFacesV

    Adding the following cell zones:
    1:  movingCellsPiston
    nV:  movingCellsPistonV

    Adding the following point zones:
    1: pistonPoints
    nV: valvePistonPointsV

*/

#   include "addPistonFacesPointZones.H"

/*
    Adding the following face zones:

    nV: curtainCylZoneV
    nV: curtainPortZoneV
    nV: cutFaceZoneV
    nV: poppetZoneV
    nV: bottomZoneV

    Adding the following point zones:

    nV: cutPointsV

*/


#   include "addValvesFacesPointZones.H"

/*

    Adding the following point zones:

    nV: valveTopPointsV
    nV: valveBottomPointsV

    Adding the following cell zones:

    nV: movingCellsTopV
    nV: movingCellsBotV

*/

#   include "addValvePistonCellZones.H"

#   include "addAttachDetachFaces.H"

    Info<< "Adding " << nPointZones << " point, "
        << nFaceZones << " face zones and "
        << nCellZones << " cell zones" << endl;


    pz.setSize(nPointZones);
    Info << "setSize pz" << endl;
    fz.setSize(nFaceZones);
    Info << "setSize fz" << endl;
    cz.setSize(nCellZones);
    Info << "setSize cz" << endl;

    addZones(pz, fz, cz);

#   include "addMeshModifiers.H"

    // Calculating the virtual positions of piston and valves

    setVirtualPositions();

    // Write mesh
    write();

    Info << "virtualPistonPosition = " << virtualPistonPosition() << endl;
    Info << "piston position = " << pistonPosition() << endl;
}


// ************************************************************************* //
