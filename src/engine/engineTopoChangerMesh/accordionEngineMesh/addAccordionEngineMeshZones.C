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

\*---------------------------------------------------------------------------*/

#include "accordionEngineMesh.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "regionSide.H"
#include "attachDetachFunctions.H"
#include "directTopoChange.H"
#include "zeroGradientTetPolyPatchFields.H"
#include "tetMotionSolver.H"

#include "fixedValueTetPolyPatchFields.H"
#include "mixedTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "zeroGradientTetPolyPatchFields.H"
#include "meshTools.H"
#include "syncTools.H"
#include "HashSet.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::accordionEngineMesh::addMeshZones()
{
    // Add the zones and mesh modifiers to operate piston motion
    Switch addBoundary(true);

    if
    (
        pointZones().size() > 0
    )
    {
        Info<< "Time = " << engTime().theta() << endl;
        Info<< "void Foam::accordionEngineMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        checkAndCalculate();

        Info << "Point zones found = " << pointZones().size() << endl;
        Info << "Face zones found = " << faceZones().size() << endl;
        Info << "Cell zones found = " << cellZones().size() << endl;

    }
    else
    {
        pointZones().setSize(0);
        faceZones().setSize(0);
        cellZones().setSize(0);
        topoChanger_.setSize(0);

        Info << "checkAndCalculate()" << endl;
        checkAndCalculate();

        Info<< "Time = " << engTime().theta() << endl
            << "Adding zones to the engine mesh" << endl;

        DynamicList<pointZone*> pz;

        DynamicList<faceZone*> fz;

        DynamicList<cellZone*> cz;

        label nPointZones = 0;
        label nFaceZones = 0;
        label nCellZones = 0;


#       include "addPistonFacesPointZonesAccordionEngineMesh.H"

#       include "addValveFaceZonesAccordionEngineMesh.H"

#       include "addOutputCellsAccordionEngineMesh.H"


        Info<< "Adding " << nPointZones << " point, "
            << nFaceZones << " face zones and "
            << nCellZones << " cell zones" << endl;

        Info << boundary().size() << endl;

        pz.setSize(nPointZones);
        Info << "setSize pz" << endl;
        fz.setSize(nFaceZones);
        Info << "setSize fz" << endl;
        cz.setSize(nCellZones);
        Info << "setSize cz" << endl;

        addZones(pz, fz, cz);
    }


    // Calculating the virtual positions of piston and valves

 //   setVirtualPositions();

    // Write mesh and modifiers
    write();
    Info << "virtualPistonPosition = " << virtualPistonPosition() << endl;
    Info << "piston position = " << pistonPosition() << endl;
}


// ************************************************************************* //
