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

#include "triSurface.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "boundBox.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::triSurface, 0);


fileName triSurface::triSurfInstance(const Time& d)
{
    fileName foamName(d.caseName() + ".ftr");

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            if (file(d.path()/ts[j].name()/typeName/foamName))
            {
                if (debug)
                {
                    Pout<< " triSurface::triSurfInstance(const Time& d)"
                        << "reading " << foamName
                        << " from " << ts[j].name()/typeName
                        << endl;
                }

                return ts[j].name();
            }
        }
    }

    if (debug)
    {
        Pout<< " triSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }
    return "constant";
}


List<labelledTri> triSurface::convertToTri
(
    const faceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, faceI)
    {
        const face& f = faces[faceI];

        if (f.size() != 3)
        {
            FatalErrorIn
            (
                "triSurface::convertToTri"
                "(const faceList&, const label)"
            )   << "Face at position " << faceI
                << " does not have three vertices:" << f
                << abort(FatalError);
        }

        labelledTri& tri = triFaces[faceI];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


List<labelledTri> triSurface::convertToTri
(
    const triFaceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, faceI)
    {
        const triFace& f = faces[faceI];

        labelledTri& tri = triFaces[faceI];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void triSurface::printTriangle
(
    Ostream& os,
    const string& pre,
    const labelledTri& f,
    const pointField& points
)
{
    os 
        << pre.c_str() << "vertex numbers:"
        << f[0] << ' ' << f[1] << ' ' << f[2] << endl
        << pre.c_str() << "vertex coords :"
        << points[f[0]] << ' ' << points[f[1]] << ' ' << points[f[2]]
        << pre.c_str() << "region        :" << f.region() << endl
        << endl;
}


string triSurface::getLineNoComment(IFstream& is)
{
    string line;
    do
    {
        is.getLine(line);
    }
    while((line.size() == 0 || line[0] == '#') && is.good());

    return line;
}


// Remove non-triangles, double triangles.
void triSurface::checkTriangles(const bool verbose)
{
    // Simple check on indices ok.
    const label maxPointI = points().size() - 1;

    forAll(*this, faceI)
    {
        const labelledTri& f = (*this)[faceI];

        if
        (
            (f[0] < 0) || (f[0] > maxPointI)
         || (f[1] < 0) || (f[1] > maxPointI)
         || (f[2] < 0) || (f[2] > maxPointI)
        )
        {
            FatalErrorIn("triSurface::checkTriangles(bool)")
                << "triangle " << f
                << " uses point indices outside point range 0.."
                << maxPointI
                << exit(FatalError);
        }
    }

    // Two phase process
    //   1. mark invalid faces
    //   2. pack
    // Done to keep numbering constant in phase 1

    // List of valid triangles
    boolList valid(size(), true);
    bool hasInvalid = false;

    const labelListList& fFaces = faceFaces();

    forAll(*this, faceI)
    {
        const labelledTri& f = (*this)[faceI];

        if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
        {
            // 'degenerate' triangle check
            valid[faceI] = false;
            hasInvalid = true;

            if (verbose)
            {
                WarningIn
                (
                    "triSurface::checkTriangles(bool verbose)"
                )   << "triangle " << faceI
                    << " does not have three unique vertices:\n";
                printTriangle(Warning, "    ", f, points());
            }
        }
        else
        {
            // duplicate triangle check
            const labelList& neighbours = fFaces[faceI];

            // Check if faceNeighbours use same points as this face.
            // Note: discards normal information - sides of baffle are merged.
            forAll(neighbours, neighbourI)
            {
                if (neighbours[neighbourI] <= faceI)
                {
                    // lower numbered faces already checked
                    continue;
                }

                const labelledTri& n = (*this)[neighbours[neighbourI]];

                if
                (
                    ((f[0] == n[0]) || (f[0] == n[1]) || (f[0] == n[2]))
                 && ((f[1] == n[0]) || (f[1] == n[1]) || (f[1] == n[2]))
                 && ((f[2] == n[0]) || (f[2] == n[1]) || (f[2] == n[2]))
                )
                {
                    valid[faceI] = false;
                    hasInvalid = true;

                    if (verbose)
                    {
                        WarningIn
                        (
                            "triSurface::checkTriangles(bool verbose)"
                        )   << "triangles share the same vertices:\n"
                            << "    face 1 :" << faceI << endl;
                        printTriangle(Warning, "    ", f, points()); 

                        Warning
                            << endl
                            << "    face 2 :"
                            << neighbours[neighbourI] << endl;
                        printTriangle(Warning, "    ", n, points()); 
                    }

                    break;
                }
            }
        }
    }

    if (hasInvalid)
    {
        // Pack
        label newFaceI = 0;
        forAll(*this, faceI)
        {
            if (valid[faceI])
            {
                const labelledTri& f = (*this)[faceI];
                (*this)[newFaceI++] = f;
            }
        }

        if (verbose)
        {
            WarningIn
            (
                "triSurface::checkTriangles(bool verbose)"
            )   << "Removing " << size() - newFaceI
                << " illegal faces." << endl;
        }
        (*this).setSize(newFaceI);

        // Topology can change because of renumbering
        clearOut();
    }
}


// Check/fix edges with more than two triangles
void triSurface::checkEdges(const bool verbose)
{
    const labelListList& eFaces = edgeFaces();

    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 0)
        {
            FatalErrorIn("triSurface::checkEdges(bool verbose)")
                << "Edge " << edgeI << " with vertices " << edges()[edgeI]
                << " has no edgeFaces"
                << exit(FatalError);
        }
        else if (myFaces.size() > 2)
        {
            WarningIn
            (
                "triSurface::checkEdges(bool verbose)"
            )   << "Edge " << edgeI << " with vertices " << edges()[edgeI]
                << " has more than 2 faces connected to it : " << myFaces
                << endl;
        }
    }
}


// Check normals and orientation
boolList triSurface::checkOrientation(const bool verbose)
{
    const labelListList& eFaces = edgeFaces();
    
    // Storage for holding status of edge. True if normal flips across this
    // edge
    boolList borderEdge(nEdges(), false);

    forAll(eFaces, edgeI)
    {
        const labelList& neighbours = eFaces[edgeI];

        if (neighbours.size() == 2)
        {
            // Two triangles, A and B. Check if edge orientation is
            // anticlockwise on both.

            const labelledTri& triA = operator[](neighbours[0]);
            const labelledTri& triB = operator[](neighbours[1]);

            // Get first shared point
            label sharedA = -1;
            label sharedB = -1;
            forAll(triA, i)
            {
                sharedB = findIndex(triB, triA[i]);
                if (sharedB != -1)
                {
                    sharedA = i;
                    break;
                }
            }

            if
            (
                sharedA == -1
             || sharedB == -1
             || triA[sharedA] != triB[sharedB]
            )
            {
                FatalErrorIn("triSurface::checkOrientation(bool)")
                    << "Problem triA:" << triA << " triB:" << triB
                    << abort(FatalError);
            }

            // Check if next point along A equals previous point along B
            // or vise versa.
            label nextA = triA[triA.fcIndex(sharedA)];
            label prevA = triA[triA.rcIndex(sharedA)];

            label nextB = triB[triB.fcIndex(sharedB)];
            label prevB = triB[triB.rcIndex(sharedB)];

            if (nextA != prevB && prevA != nextB)
            {
                borderEdge[edgeI] = true;
                if (verbose)
                {
                    WarningIn("triSurface::checkOrientation(bool)")
                        << "Triangle orientation incorrect." << endl
                        << "edge neighbours:" << neighbours << endl
                        << "triangle " << neighbours[0]
                        << " points " << triA << endl
                        << "triangle " << neighbours[1]
                        << " points " << triB << endl;
                }
            }
        }
        else if (neighbours.size() != 1)
        {
            if (verbose)
            {
                const edge& e = edges()[edgeI];
                WarningIn("triSurface::checkOrientation(bool)")
                    << "Wrong number of edge neighbours." << endl
                    << "Edge:" << e
                    << "with points:" << localPoints()[e.start()]
                    << ' ' << localPoints()[e.end()]
                    << " has neighbours:" << neighbours << endl;
            }
            borderEdge[edgeI] = true;
        }
    }

    return borderEdge;
}


// Read triangles, points from Istream
bool triSurface::read(Istream& is)
{
    is  >> patches_ >> const_cast<pointField&>(points())
        >> static_cast<List<labelledTri>&>(*this);

    return true;
}


// Read from file in given format
bool triSurface::read(const fileName& name, const word& ext, const bool check)
{
    if (check && !exists(name))
    {
        FatalErrorIn
        (
            "triSurface::read(const fileName&, const word&, const bool)"
        )   << "Cannnot read " << name << exit(FatalError);
    }

    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();

        // Do not check for existence. Let IFstream do the unzipping.
        return read(unzipName, unzipName.ext(), false);
    }
    else if (ext == "ftr")
    {
        return read(IFstream(name)());
    }
    else if (ext == "stl")
    {
        return readSTL(name);
    }
    else if (ext == "stlb")
    {
        return readSTL(name);
    }
    else if (ext == "gts")
    {
        return readGTS(name);
    }
    else if (ext == "obj")
    {
        return readOBJ(name);
    }
    else if (ext == "off")
    {
        return readOFF(name);
    }
    else if (ext == "tri")
    {
        return readTRI(name);
    }
    else if (ext == "ac")
    {
        return readAC(name);
    }
    else if (ext == "nas")
    {
        return readNAS(name);
    }
    else
    {
        FatalErrorIn
        (
            "triSurface::read(const fileName&, const word&)"
        )   << "unknown file extension " << ext
            << ". Supported extensions are '.ftr', '.stl', '.stlb', '.gts'"
            << ", '.obj', '.ac', '.off', '.nas' and '.tri'"
            << exit(FatalError);

        return false;
    }
}


// Write to file in given format
void triSurface::write
(
    const fileName& name,
    const word& ext,
    const bool sort
) const
{
    if (ext == "ftr")
    {
        return write(OFstream(name)());
    }
    else if (ext == "stl")
    {
        return writeSTLASCII(OFstream(name)());
    }
    else if (ext == "stlb")
    {
        ofstream outFile(name.c_str(), std::ios::binary);

        writeSTLBINARY(outFile);
    }
    else if (ext == "gts")
    {
        return writeGTS(sort, OFstream(name)());
    }
    else if (ext == "obj")
    {
        writeOBJ(sort, OFstream(name)());
    }
    else if (ext == "off")
    {
        writeOFF(sort, OFstream(name)());
    }
    else if (ext == "vtk")
    {
        writeVTK(sort, OFstream(name)());
    }
    else if (ext == "tri")
    {
        writeTRI(sort, OFstream(name)());
    }
    else if (ext == "dx")
    {
        writeDX(sort, OFstream(name)());
    }
    else if (ext == "ac")
    {
        writeAC(OFstream(name)());
    }
    else if (ext == "smesh")
    {
        writeSMESH(sort, OFstream(name)());
    }
    else
    {
        FatalErrorIn
        (
            "triSurface::write(const fileName&, const word&, const bool)"
        )   << "unknown file extension " << ext
            << ". Supported extensions are '.ftr', '.stl', '.stlb', "
            << "'.gts', '.obj', '.vtk'"
            << ", '.off', '.dx', '.smesh', '.ac' and '.tri'"
            << exit(FatalError);
    }
}


// Returns patch info. Sets faceMap to the indexing according to patch
// numbers. Patch numbers start at 0.
surfacePatchList triSurface::calcPatches(labelList& faceMap) const
{
    // Sort according to region numbers of labelledTri
    SortableList<label> sortedRegion(size());

    forAll(sortedRegion, faceI)
    {
        sortedRegion[faceI] = operator[](faceI).region();
    }
    sortedRegion.stableSort();

    faceMap = sortedRegion.indices();

    // Compact regions

    // Get last region
    label maxRegion = 0;

    if (faceMap.size() > 0)
    {
        maxRegion = operator[](faceMap[faceMap.size() - 1]).region();
    }

    // Get new region list
    surfacePatchList newPatches(maxRegion + 1);

    // Fill patch sizes
    forAll(*this, faceI)
    {
        label region = operator[](faceI).region();

        newPatches[region].size()++;
    }


    // Fill rest of patch info

    label startFaceI = 0;
    forAll(newPatches, newPatchI)
    {
        surfacePatch& newPatch = newPatches[newPatchI];

        label oldPatchI = newPatchI;

        // start of patch
        newPatch.start() = startFaceI;


        // Take over any information from existing patches
        if ((oldPatchI < patches_.size()) && (patches_[oldPatchI].name() != ""))
        {
            newPatch.name() = patches_[oldPatchI].name();
        }
        else
        {
            newPatch.name() = word("patch") + name(newPatchI);
        }

        if
        (
            (oldPatchI < patches_.size())
         && (patches_[oldPatchI].geometricType() != "")
        )
        {
            newPatch.geometricType() = patches_[oldPatchI].geometricType();
        }
        else
        {
            newPatch.geometricType() = "empty";
        }

        startFaceI += newPatch.size();
    }

    return newPatches;
}


void triSurface::setDefaultPatches()
{
    labelList faceMap;

    // Get names, types and sizes
    surfacePatchList newPatches(calcPatches(faceMap));

    // Take over names and types (but not size)
    patches_.setSize(newPatches.size());

    forAll(newPatches, patchI)
    {
        patches_[patchI].name() = newPatches[patchI].name();
        patches_[patchI].geometricType() = newPatches[patchI].geometricType();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurface::triSurface()
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>
    (
        List<labelledTri>(0),
        pointField(0)
    ),
    patches_(0),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{}



triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const pointField& points
)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>(triangles, points),
    patches_(patches),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{}

triSurface::triSurface
(
    List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    pointField& points,
    const bool reUse
)
:
    PrimitivePatch<labelledTri, ::Foam::List ,pointField>(triangles, points, reUse),
    patches_(patches),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{}

triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const pointField& points
)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>(triangles, points),
    patches_(),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{
    setDefaultPatches();
}


triSurface::triSurface
(
    const triFaceList& triangles,
    const pointField& points
)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>
    (
        convertToTri(triangles, 0),
        points
    ),
    patches_(),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{
    setDefaultPatches();
}


triSurface::triSurface(const fileName& name)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>
    (
        List<labelledTri>(0),
        pointField(0)
    ),
    patches_(),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{
    word ext = name.ext();

    read(name, ext);

    setDefaultPatches();
}


triSurface::triSurface(Istream& is)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>
    (
        List<labelledTri>(0),
        pointField(0)
    ),
    patches_(),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{
    read(is);

    setDefaultPatches();
}


triSurface::triSurface(const Time& d)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>
    (
        List<labelledTri>(0),
        pointField(0)
    ),
    patches_(),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{
    fileName foamFile(d.caseName() + ".ftr");

    fileName foamPath(d.path()/triSurfInstance(d)/typeName/foamFile);

    IFstream foamStream(foamPath);

    read(foamStream);

    setDefaultPatches();
}


triSurface::triSurface(const triSurface& ts)
:
    PrimitivePatch<labelledTri, ::Foam::List, pointField>(ts, ts.points()),
    patches_(ts.patches()),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triSurface::~triSurface()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::clearTopology()
{
    PrimitivePatch<labelledTri, ::Foam::List, pointField>::clearTopology();
    deleteDemandDrivenData(sortedEdgeFacesPtr_);
    deleteDemandDrivenData(edgeOwnerPtr_);
}


void triSurface::clearPatchMeshAddr()
{
    PrimitivePatch<labelledTri, ::Foam::List, pointField>::clearPatchMeshAddr();
}


void triSurface::clearOut()
{
    PrimitivePatch<labelledTri, ::Foam::List, pointField>::clearOut();

    clearTopology();
    clearPatchMeshAddr();
}


const labelListList& triSurface::sortedEdgeFaces() const
{
    if (!sortedEdgeFacesPtr_)
    {
        calcSortedEdgeFaces();
    }

    return *sortedEdgeFacesPtr_;
}


const labelList& triSurface::edgeOwner() const
{
    if (!edgeOwnerPtr_)
    {
        calcEdgeOwner();
    }

    return *edgeOwnerPtr_;
}


//- Move points
void triSurface::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    deleteDemandDrivenData(sortedEdgeFacesPtr_);

    // Adapt for new point position
    PrimitivePatch<labelledTri, ::Foam::List, pointField>::movePoints(newPoints);

    // Copy new points
    const_cast<pointField&>(points()) = newPoints;
}

// scale points
void triSurface::scalePoints(const scalar& scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        // Remove all geometry dependent data
        clearTopology();

        // Adapt for new point position
        PrimitivePatch<labelledTri, ::Foam::List, pointField>::movePoints(pointField());

        const_cast<pointField&>(points()) *= scaleFactor;
    }
}

// Remove non-triangles, double triangles.
void triSurface::cleanup(const bool verbose)
{
    // Merge points (already done for STL, TRI)
    stitchTriangles(pointField(points()), SMALL, verbose);

    // Merging points might have changed geometric factors
    clearOut();

    checkTriangles(verbose);

    checkEdges(verbose);
}


// Finds area, starting at faceI, delimited by borderEdge. Marks all visited
// faces (from face-edge-face walk) with currentZone.
void triSurface::markZone
(
    const boolList& borderEdge,
    const label faceI,
    const label currentZone,
    labelList& faceZone
) const
{
    // List of faces whose faceZone has been set.
    labelList changedFaces(1, faceI);

    while(true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        forAll(changedFaces, i)
        {
            label faceI = changedFaces[i];

            const labelList& fEdges = faceEdges()[faceI];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                if (!borderEdge[edgeI])
                {
                    const labelList& eFaces = edgeFaces()[edgeI];

                    forAll(eFaces, j)
                    {
                        label nbrFaceI = eFaces[j];

                        if (faceZone[nbrFaceI] == -1)
                        {
                            faceZone[nbrFaceI] = currentZone;
                            newChangedFaces.append(nbrFaceI);
                        }
                        else if (faceZone[nbrFaceI] != currentZone)
                        {
                            FatalErrorIn
                            (
                                "triSurface::markZone(const boolList&,"
                                "const label, const label, labelList&) const"
                            )
                                << "Zones " << faceZone[nbrFaceI]
                                << " at face " << nbrFaceI
                                << " connects to zone " << currentZone
                                << " at face " << faceI
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.size() == 0)
        {
            break;
        }

        changedFaces.transfer(newChangedFaces.shrink());
        newChangedFaces.clear();
    }
}


// Finds areas delimited by borderEdge (or 'real' edges).
// Fills faceZone accordingly
label triSurface::markZones
(
    const boolList& borderEdge,
    labelList& faceZone
) const
{
    faceZone.setSize(size());
    faceZone = -1;

    if (borderEdge.size() != nEdges())
    {
        FatalErrorIn
        (
            "triSurface::markZones"
            "(const boolList&, labelList&)"
        )
            << "borderEdge boolList not same size as number of edges" << endl
            << "borderEdge:" << borderEdge.size() << endl
            << "nEdges    :" << nEdges()
            << exit(FatalError);
    }

    label zoneI = 0;

    label startFaceI = 0;

    for(;;zoneI++)
    {
        // Find first non-coloured face
        for (; startFaceI < size(); startFaceI++)
        {
            if (faceZone[startFaceI] == -1)
            {
                break;
            }
        }

        if (startFaceI >= size())
        {
            break;
        }

        faceZone[startFaceI] = zoneI;

        markZone(borderEdge, startFaceI, zoneI, faceZone);
    }

    return zoneI;
}


void triSurface::subsetMeshMap
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const List<labelledTri>& locFaces = localFaces();

    label faceI = 0;
    label pointI = 0;

    faceMap.setSize(locFaces.size());

    pointMap.setSize(nPoints());

    boolList pointHad(nPoints(), false);

    forAll(include, oldFacei)
    {
        if (include[oldFacei])
        {
            // Store new faces compact
            faceMap[faceI++] = oldFacei;

            // Renumber labels for triangle
            const labelledTri& tri = locFaces[oldFacei];

            label a = tri[0];
            if (!pointHad[a])
            {
                pointHad[a] = true;
                pointMap[pointI++] = a;                
            }

            label b = tri[1];
            if (!pointHad[b])
            {
                pointHad[b] = true;
                pointMap[pointI++] = b;
            }

            label c = tri[2];
            if (!pointHad[c])
            {
                pointHad[c] = true;
                pointMap[pointI++] = c;
            }
        }
    }

    // Trim
    faceMap.setSize(faceI);

    pointMap.setSize(pointI);
}


triSurface triSurface::subsetMesh
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = localPoints();
    const List<labelledTri>& locFaces = localFaces();

    // Fill pointMap, faceMap
    subsetMeshMap(include, pointMap, faceMap);
    

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointi)
    {
        newPoints[pointi] = locPoints[pointMap[pointi]];
        oldToNew[pointMap[pointi]] = pointi;
    }

    // Renumber triangle node labels and compact
    List<labelledTri> newTriangles(faceMap.size());

    forAll(faceMap, facei)
    {
        // Get old vertex labels
        const labelledTri& tri = locFaces[faceMap[facei]];

        newTriangles[facei][0] = oldToNew[tri[0]];
        newTriangles[facei][1] = oldToNew[tri[1]];
        newTriangles[facei][2] = oldToNew[tri[2]];
        newTriangles[facei].region() = tri.region();
    }

    // Construct subsurface
    triSurface subSurface(newTriangles, patches(), newPoints);

    return subSurface;
}


void triSurface::write(const fileName& name, const bool sortByRegion) const
{
    write(name, name.ext(), sortByRegion);
}


void triSurface::write(Ostream& os) const
{
    os  << patches() << endl;

    //Note: Write with global point numbering
    os  << points() << nl
        << static_cast<const List<labelledTri>&>(*this) << endl;

    // Check state of Ostream
    os.check("triSurface::write(Ostream&)");
}


void triSurface::write(const Time& d) const
{
    fileName foamFile(d.caseName() + ".ftr");

    fileName foamPath(d.path()/triSurfInstance(d)/typeName/foamFile);

    OFstream foamStream(foamPath);

    write(foamStream);
}


void triSurface::writeStats(Ostream& os) const
{
    // Calculate bounding box without any additional addressing
    // Copy of treeBoundBox code. Cannot use meshTools from triSurface...
    boundBox bb
    (
        point(VGREAT, VGREAT, VGREAT),
        point(-VGREAT, -VGREAT, -VGREAT)
    );
    forAll(*this, triI)
    {
        const labelledTri& f = operator[](triI);

        forAll(f, fp)
        {
            const point& pt = points()[f[fp]];
            bb.min() = ::Foam::min(bb.min(), pt);
            bb.max() = ::Foam::max(bb.max(), pt);
        }
    }

    // Unfortunately nPoints constructs meshPoints() ...

    os  << "Triangles    : " << size() << endl
        //<< "Edges        : " << nEdges() << endl
        << "Vertices     : " << nPoints() << endl
        << "Bounding Box : " << bb << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void triSurface::operator=(const triSurface& ts)
{
    List<labelledTri>::operator=(ts);
    clearOut();
    const_cast<pointField&>(points()) = ts.points();
    patches_ = ts.patches();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const triSurface& sm)
{
    sm.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
