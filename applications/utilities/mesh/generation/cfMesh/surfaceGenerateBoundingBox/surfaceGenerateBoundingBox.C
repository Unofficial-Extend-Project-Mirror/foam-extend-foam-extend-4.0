/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description
    Finds feature edges and corners of a triangulated surface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "boundBox.H"
#include "OFstream.H"

#include <cstdlib>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("x-neg");
    argList::validArgs.append("x-pos");
    argList::validArgs.append("y-neg");
    argList::validArgs.append("y-pos");
    argList::validArgs.append("z-neg");
    argList::validArgs.append("z-pos");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if (outFileName == inFileName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite the input file."
            << exit(FatalError);
    }

    triSurf origSurface(inFileName);
    triSurfModifier sMod(origSurface);
    pointField& points = sMod.pointsAccess();

    const boundBox bb(points);

    vector negOffset, posOffset;
    for(label i=3;i<9;++i)
    {
        std::stringstream ss;
        ss << args.args()[i];

        scalar s;
        ss >> s;

        if( i % 2 )
        {
            negOffset[(i-3)/2] = s;
        }
        else
        {
            posOffset[(i-3)/2] = s;
        }
    }

    Info << "Neg offset " << negOffset << endl;
    Info << "Pos offset " << posOffset << endl;

    const boundBox newBB(bb.min()-negOffset, bb.max()+posOffset);
    Info << "Surface bounding box " << bb << endl;
    Info << "Generated bounding box " << newBB << endl;

    //- generate bounding box points
    const label nPoints = points.size();
    points.setSize(nPoints + 8);

    points[nPoints] = newBB.min();
    points[nPoints+1] =
        point(newBB.max().x(), newBB.min().y(), newBB.min().z());
    points[nPoints+2] =
        point(newBB.min().x(), newBB.max().y(), newBB.min().z());
    points[nPoints+3] =
        point(newBB.max().x(), newBB.max().y(), newBB.min().z());
    points[nPoints+4] =
        point(newBB.min().x(), newBB.min().y(), newBB.max().z());
    points[nPoints+5] =
        point(newBB.max().x(), newBB.min().y(), newBB.max().z());
    points[nPoints+6] =
        point(newBB.min().x(), newBB.max().y(), newBB.max().z());
    points[nPoints+7] = newBB.max();

    //- generate bounding bound triangles
    const label nTriangles = origSurface.size();
    LongList<labelledTri>& newTriangles = sMod.facetsAccess();
    newTriangles.setSize(nTriangles+12);

    //- create patches
    geometricSurfacePatchList& newPatches = sMod.patchesAccess();
    const label nPatches = origSurface.patches().size();
    newPatches.setSize(nPatches+6);

    newPatches[nPatches].name() = "xMin";
    newPatches[nPatches+1].name() = "xMax";
    newPatches[nPatches+2].name() = "yMin";
    newPatches[nPatches+3].name() = "yMax";
    newPatches[nPatches+4].name() = "zMin";
    newPatches[nPatches+5].name() = "zMax";

    //- negative x direction
    newTriangles[nTriangles] =
        labelledTri(nPoints, nPoints+6, nPoints+2, nPatches);
    newTriangles[nTriangles+1] =
        labelledTri(nPoints, nPoints+4, nPoints+6, nPatches);
    //- positive x direction
    newTriangles[nTriangles+2] =
        labelledTri(nPoints+1, nPoints+3, nPoints+7, nPatches+1);
    newTriangles[nTriangles+3] =
        labelledTri(nPoints+1, nPoints+7, nPoints+5, nPatches+1);
    //- negative y direction
    newTriangles[nTriangles+4] =
        labelledTri(nPoints, nPoints+1, nPoints+5, nPatches+2);
    newTriangles[nTriangles+5] =
        labelledTri(nPoints, nPoints+5, nPoints+4, nPatches+2);
    //- positive y direction
    newTriangles[nTriangles+6] =
        labelledTri(nPoints+2, nPoints+7, nPoints+3, nPatches+3);
    newTriangles[nTriangles+7] =
        labelledTri(nPoints+2, nPoints+6, nPoints+7, nPatches+3);
    //- negative z direction
    newTriangles[nTriangles+8] =
        labelledTri(nPoints, nPoints+2, nPoints+3, nPatches+4);
    newTriangles[nTriangles+9] =
        labelledTri(nPoints, nPoints+3, nPoints+1, nPatches+4);
    //- positive z direction
    newTriangles[nTriangles+10] =
        labelledTri(nPoints+4, nPoints+7, nPoints+6, nPatches+5);
    newTriangles[nTriangles+11] =
        labelledTri(nPoints+4, nPoints+5, nPoints+7, nPatches+5);

    //- write the surface
    origSurface.writeSurface(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
