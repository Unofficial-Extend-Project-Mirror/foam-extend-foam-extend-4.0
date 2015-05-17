/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

*/
// ------------------------------------------------------------------------- //

#include "pistonLayer.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "pointSet.H"
#include "faceSet.H"
#include "SortableList.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pistonLayer::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston motion

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void pistonLayer::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void pistonLayer::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        Info << "topoChanger.size() = " << topoChanger_.size() << endl;

        checkAndCalculate();

        setVirtualPistonPosition();
        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
        topoChanger_.write();

        return;
    }

    checkAndCalculate();

    Info<< "Time = " << engTime().theta() << endl
        << "Adding zones to the engine mesh" << endl;

    //fz = 1: faces where layer are added/removed
    //pz = 2: points below the virtual piston faces and head points

    DynamicList<pointZone*> pz(2);
    DynamicList<faceZone*> fz(1);
    DynamicList<cellZone*> cz(0);

    label nPointZones = 0;
    label nFaceZones = 0;

#   include "addPistonLayerFaces.H"

    {

        pointSet headPointSet(*this, headPointsSetName_);

        Info << "Number of head points = " << headPointSet.size() << endl;
        pz[nPointZones] =
            new pointZone
            (
                "headPoints",
                headPointSet.toc(),
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

#   include "addPistonLayerAdditionRemovalMeshModifier.H"

    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();

    write();

    // Calculating the virtual piston position
    setVirtualPistonPosition();

    Info << "virtualPistonPosition = " << virtualPistonPosition() << endl;
    Info << "piston position = " << pistonPosition() << endl;

}


// ************************************************************************* //
