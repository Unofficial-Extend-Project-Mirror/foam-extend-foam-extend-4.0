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
    Check a Finite Area mesh

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFaMesh.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // General mesh statistics
    Info<< "Number of points: " << aMesh.nPoints() << nl
        << "Number of internal edges: " << aMesh.nInternalEdges() << nl
        << "Number of edges: " << aMesh.nEdges() << nl
        << "Number of faces: " << aMesh.nFaces() << nl
        << endl;

    // Check geometry
    Info<< "Face area:  min = " << min(aMesh.S().field())
        << " max = "  << max(aMesh.S().field()) << nl
        << "Internal edge length: min = "
        << min(aMesh.magLe().internalField()) << nl
        << " max = "  << max(aMesh.magLe().internalField()) << nl
        << "Edge length: min = "
        << min(aMesh.magLe()).value() << nl
        << " max = "  << max(aMesh.magLe()).value() << nl
        << "Face area normals:  min = " << min(aMesh.faceAreaNormals().field())
        << " max = "  << max(aMesh.faceAreaNormals().field()) << nl
        << endl;


    return(0);
}


// ************************************************************************* //
