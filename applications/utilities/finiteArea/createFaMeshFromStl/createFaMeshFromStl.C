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

Application
    Create a Finite Area mesh from an STL file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "triSurface.H"
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
