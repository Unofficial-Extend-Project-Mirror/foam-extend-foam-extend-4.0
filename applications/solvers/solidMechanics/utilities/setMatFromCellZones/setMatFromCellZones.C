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
    setMatFromCellZones

Description
    Reads in the mesh cellZones and then create a material volScalarField
    based on the cellZones.
    This is works well for a multi-zone mesh created in snappyHexMesh.

Author
   Philip Cardiff

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMesh.H"

    Info<< "Reading mesh cellZones" << endl;

    const cellZoneMesh& cellZones = mesh.cellZones();

    forAll(cellZones, zonei)
      {
	Info << "\tCell zone " << cellZones[zonei].name()
	     << " with " << cellZones[zonei].size() << " cells"
	     << endl;
      }

    Info << "\nCreating materials field\n" << endl;

    volScalarField materials
    (
        IOobject
        (
            "materials",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellZonesSize", dimless, cellZones.size()),
        zeroGradientFvPatchScalarField::typeName
    );
    //- the default value of materials is set to be the size of cellZones
    //- so cellZone[0] will be set to 0, cellZone[1] will be set to 1,
    //- and so on and then any cell not in a cellZone will be left as
    //- cellZones.size()
    scalarField& materialsI = materials.internalField();
    
    forAll(cellZones, zonei)
      {
	labelList zoneCells = cellZones[zonei];

	forAll(zoneCells, celli)
	  {
	    const label& ci = zoneCells[celli];
	    materialsI[ci] = zonei;
	  }
      }

    materials.correctBoundaryConditions();

    //- check all cells have been set
    if(max(materialsI) == cellZones.size())
      {
	Warning << "There are cells which are not in a cellZone"
		", their material field is set to "
		<< cellZones.size() << endl;
      }

    Info << "\nWriting materials field\n" << endl;
    materials.write();

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
