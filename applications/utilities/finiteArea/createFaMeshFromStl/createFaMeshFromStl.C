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
    Create a Finite Area mesh from an STL file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"
#include "Time.H"
#include "polyMesh.H"
#include "faCFD.H"
#include "IFstream.H"
#include "graph.H"
#include "Tuple2.H"
#include "matchPoints.H"
#include "standAlonePatch.H"

#include "fixedValueFaPatchFields.H"
#include "zeroGradientFaPatchFields.H"
#include "inletOutletFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("STL mesh file");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    fileName prefix(args.additionalArgs()[0]);

#   include "createTime.H"
#   include "makePolyMesh.H"
#   include "makeFaMesh.H"

    Info<< " done." << nl << endl;

    Info<< nl << "Write mesh vtk files... ";
    standAlonePatch::writeVTK
    (
        runTime.caseName() + "Mesh",
        aMesh.patch().localFaces(),
        aMesh.patch().localPoints()
    );

    standAlonePatch::writeVTKNormals
    (
        runTime.caseName() + "Normals",
        aMesh.patch().localFaces(),
        aMesh.patch().localPoints()
    );

    Info<< " done." << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
