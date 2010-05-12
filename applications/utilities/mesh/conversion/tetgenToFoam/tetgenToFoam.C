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

Description
    Reads .ele and .node and .face files as written by tetgen.

    Make sure to use add boundary attributes to the smesh file
    (5 fifth column in the element section)
    and run tetgen with -f option.

    Sample smesh file:

        # cube.smesh -- A 10x10x10 cube
        8 3
        1	0 0 0
        2	0 10 0
        3	10 10 0
        4	10 0 0
        5	0 0 10
        6	0 10 10
        7	10 10 10
        8	10 0 10
        6 1                 # 1 for boundary info present
        4	1 2 3 4 11  # region number 11
        4	5 6 7 8 21  # region number 21
        4	1 2 6 5 3
        4	4 3 7 8 43
        4	1 5 8 4 5
        4	2 6 7 3 65
        0
        0

NOTE:
- for some reason boundary faces point inwards. I just reverse them
always. Might use some geometric check instead.
- marked faces might not actually be boundary faces of mesh. This is not handled
and you'll have to run without face file (-noFaceFile option)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "cellModeller.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("file prefix");
    argList::validOptions.insert("noFaceFile", "");
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"


    bool readFaceFile = !args.options().found("noFaceFile");
    bool overwrite = args.options().found("overwrite");

    fileName prefix(args.additionalArgs()[0]);

    fileName nodeFile(prefix + ".node");
    fileName eleFile(prefix + ".ele");
    fileName faceFile(prefix + ".face");

    if (!readFaceFile)
    {
        Info<< "Files:" << endl
            << "    nodes : " << nodeFile << endl
            << "    elems : " << eleFile << endl
            << endl;
    }
    else
    {
        Info<< "Files:" << endl
            << "    nodes : " << nodeFile << endl
            << "    elems : " << eleFile << endl
            << "    faces : " << faceFile << endl
            << endl;

        Info<< "Reading .face file for boundary information" << nl << endl;
    }

    if (!exists(nodeFile) || !exists(eleFile))
    {
        FatalErrorIn(args.executable())
            << "Cannot read " << nodeFile << " or " << eleFile
            << exit(FatalError);
    }

    if (readFaceFile && !exists(faceFile))
    {
        FatalErrorIn(args.executable())
            << "Cannot read " << faceFile << endl
            << "Did you run tetgen with -f option?" << endl
            << "If you don't want to read the .face file and thus not have"
            << " patches please\nrerun with the -noFaceFile option"
            << exit(FatalError);
    }


    IFstream nodeStream(nodeFile);

    //
    // Read nodes.
    //

    // Read header.
    string line;

    do
    {
        nodeStream.getLine(line);
    }
    while((line.size() > 0) && (line[0] == '#'));

    IStringStream nodeLine(line);

    label nNodes, nDims, nNodeAttr;
    bool hasRegion;

    nodeLine >> nNodes >> nDims >> nNodeAttr >> hasRegion;


    Info<< "Read .node header:" << endl
        << "    nodes     : " << nNodes << endl
        << "    nDims     : " << nDims << endl
        << "    nAttr     : " << nNodeAttr << endl
        << "    hasRegion : " << hasRegion << endl
        << endl;

    //
    // read points
    //

    pointField points(nNodes);
    Map<label> nodeToPoint(nNodes);

    {
        labelList pointIndex(nNodes);

        label pointI = 0;

        while (nodeStream.good())
        {
            nodeStream.getLine(line);

            if ((line.size() > 0) && (line[0] != '#'))
            {
                IStringStream nodeLine(line);

                label nodeI;
                scalar x, y, z;
                label dummy;

                nodeLine >> nodeI >> x >> y >> z;

                for (label i = 0; i < nNodeAttr; i++)
                {
                    nodeLine >> dummy;
                }

                if (hasRegion)
                {
                    nodeLine >> dummy;
                }

                // Store point and node number.
                points[pointI] = point(x, y, z);
                nodeToPoint.insert(nodeI, pointI);
                pointI++;
            }
        }
        if (pointI != nNodes)
        {
            FatalIOErrorIn(args.executable().c_str(), nodeStream)
                << "Only " << pointI << " nodes present instead of " << nNodes
                << " from header." << exit(FatalIOError);
        }
    }

    //
    // read elements
    //

    IFstream eleStream(eleFile);

    do
    {
        eleStream.getLine(line);
    }
    while((line.size() > 0) && (line[0] == '#'));

    IStringStream eleLine(line);

    label nTets, nPtsPerTet, nElemAttr;

    eleLine >> nTets >> nPtsPerTet >> nElemAttr;


    Info<< "Read .ele header:" << endl
        << "    tets         : " << nTets << endl
        << "    pointsPerTet : " << nPtsPerTet << endl
        << "    nAttr        : " << nElemAttr << endl
        << endl;

    if (nPtsPerTet != 4)
    {
        FatalIOErrorIn(args.executable().c_str(), eleStream)
            << "Cannot handle tets with "
            << nPtsPerTet << " points per tetrahedron in .ele file" << endl
            << "Can only handle tetrahedra with four points"
            << exit(FatalIOError);
    }

    if (nElemAttr != 0)
    {
        WarningIn(args.executable())
            << "Element attributes (third elemenent in .ele header)"
            << " not used" << endl;
    }



    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelList tetPoints(4);

    cellShapeList cells(nTets);
    label cellI = 0;

    while (eleStream.good())
    {
        eleStream.getLine(line);

        if ((line.size() > 0) && (line[0] != '#'))
        {
            IStringStream eleLine(line);

            label elemI;
            eleLine >> elemI;

            for (label i = 0; i < 4; i++)
            {
                label nodeI;
                eleLine >> nodeI;
                tetPoints[i] = nodeToPoint[nodeI];
            }

            cells[cellI++] = cellShape(tet, tetPoints);

            // Skip attributes
            for (label i = 0; i < nElemAttr; i++)
            {
                label dummy;

                eleLine >> dummy;
            }
        }
    }


    label nPatches = 0;

    // List of Foam vertices per boundary face
    faceList boundaryFaces;

    // For each boundary faces the Foam patchID
    labelList boundaryPatch;

    if (readFaceFile)
    {
        //
        // read boundary faces
        //

        IFstream faceStream(faceFile);

        do
        {
            faceStream.getLine(line);
        }
        while((line.size() > 0) && (line[0] == '#'));

        IStringStream faceLine(line);

        label nFaces, nFaceAttr;

        faceLine >> nFaces >> nFaceAttr;


        Info<< "Read .face header:" << endl
            << "    faces : " << nFaces << endl
            << "    nAttr : " << nFaceAttr << endl
            << endl;


        if (nFaceAttr != 1)
        {
            FatalIOErrorIn(args.executable().c_str(), faceStream)
                << "Expect boundary markers to be"
                << " present in .face file." << endl
                << "This is the second number in the header which is now:"
                << nFaceAttr << exit(FatalIOError);
        }

        // List of Foam vertices per boundary face
        boundaryFaces.setSize(nFaces);

        // For each boundary faces the Foam patchID
        boundaryPatch.setSize(nFaces);
        boundaryPatch = -1;

        label faceI = 0;

        // Region to patch conversion
        Map<label> regionToPatch;

        face f(3);

        while (faceStream.good())
        {
            faceStream.getLine(line);

            if ((line.size() > 0) && (line[0] != '#'))
            {
                IStringStream faceLine(line);

                label tetGenFaceI, dummy, region;

                faceLine >> tetGenFaceI;

                // Read face and reverse orientation (Foam needs outwards
                // pointing)
                for (label i = 0; i < 3; i++)
                {
                    label nodeI;
                    faceLine >> nodeI;
                    f[2-i] = nodeToPoint[nodeI];
                }

                boundaryFaces[faceI] = f;

                if (nFaceAttr > 0)
                {
                    // First attribute is the region number
                    faceLine >> region;


                    // Get Foam patchID and update region->patch table.
                    label patchI = 0;

                    Map<label>::iterator patchFind = regionToPatch.find(region);

                    if (patchFind == regionToPatch.end())
                    {
                        patchI = nPatches;

                        Info<< "Mapping tetgen region " << region
                            << " to Foam patch "
                            << patchI << endl;

                        regionToPatch.insert(region, nPatches++);
                    }
                    else
                    {
                        patchI = patchFind();
                    }

                    boundaryPatch[faceI] = patchI;

                    // Skip remaining attributes
                    for (label i = 1; i < nFaceAttr; i++)
                    {
                        faceLine >> dummy;
                    }
                }

                faceI++;
            }
        }

        // Print region to patch mapping
        Info<< "Regions:" << endl;

        for
        (
            Map<label>::const_iterator iter = regionToPatch.begin();
            iter != regionToPatch.end();
            ++iter
        )
        {
            Info<< "    region:" << iter.key() << '\t' << "patch:"
                << iter() << endl;
        }
        Info<< endl;
    }


    // Storage for boundary faces

    faceListList patchFaces(nPatches);

    wordList patchNames(nPatches);

    forAll(patchNames, patchI)
    {
        patchNames[patchI] = word("patch") + name(patchI);
    }

    wordList patchTypes(nPatches, polyPatch::typeName);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = polyPatch::typeName;
    wordList patchPhysicalTypes(nPatches, polyPatch::typeName);


    if (readFaceFile)
    {
        // Sort boundaryFaces by patch using boundaryPatch.
        List<DynamicList<face> > allPatchFaces(nPatches);

        forAll(boundaryPatch, faceI)
        {
            label patchI = boundaryPatch[faceI];

            allPatchFaces[patchI].append(boundaryFaces[faceI]);
        }

        Info<< "Patch sizes:" << endl;

        forAll(allPatchFaces, patchI)
        {
            Info<< "    " << patchNames[patchI] << " : "
                << allPatchFaces[patchI].size() << endl;

            patchFaces[patchI].transfer(allPatchFaces[patchI].shrink());
        }

        Info<< endl;
    }

    if (!overwrite)
    {
        runTime++;
    }

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        points,
        cells,
        patchFaces,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    Info<< "Writing mesh to " << runTime.constant() << endl << endl;

    mesh.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
