/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Class
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "foamTime.H"
#include "meshOps.H"
#include "ListOps.H"
#include "Pstream.H"
#include "triFace.H"
#include "IOmanip.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

namespace meshOps
{

// Utility method to build a hull of cells
// connected to the edge (for 2D simplical meshes)
inline void constructPrismHull
(
    const label eIndex,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    const UList<label>& neighbour,
    const UList<labelList>& edgeFaces,
    labelList& hullTriFaces,
    labelList& hullCells
)
{
    labelHashSet cellSet, triFaceSet;

    // Obtain references
    const labelList& eFaces = edgeFaces[eIndex];

    // Loop through edgeFaces and add cells
    forAll(eFaces, faceI)
    {
        label c0 = owner[eFaces[faceI]];
        label c1 = neighbour[eFaces[faceI]];

        if (!cellSet.found(c0))
        {
            // Add this cell
            cellSet.insert(c0);

            // Find associated triFaces and add them too
            const cell& cC = cells[c0];

            forAll(cC, faceJ)
            {
                const face& cF = faces[cC[faceJ]];

                if ((cF.size() == 3) && !(triFaceSet.found(cC[faceJ])))
                {
                    triFaceSet.insert(cC[faceJ]);
                }
            }
        }

        if (!cellSet.found(c1) && (c1 != -1))
        {
            // Add this cell
            cellSet.insert(c1);

            // Find associated triFaces and add them too
            const cell& cC = cells[c1];

            forAll(cC, faceJ)
            {
                const face& cF = faces[cC[faceJ]];

                if ((cF.size() == 3) && !(triFaceSet.found(cC[faceJ])))
                {
                    triFaceSet.insert(cC[faceJ]);
                }
            }
        }
    }

    // Obtain lists from hashSets
    hullCells = cellSet.toc();
    hullTriFaces = triFaceSet.toc();
}


// Utility method to build a hull of cells (and faces)
// around an edge (for 3D simplical meshes)
inline void constructHull
(
    const label eIndex,
    const UList<face>& faces,
    const UList<edge>& edges,
    const UList<cell>& cells,
    const UList<label>& owner,
    const UList<label>& neighbour,
    const UList<labelList>& faceEdges,
    const UList<labelList>& edgeFaces,
    const labelList& vertexHull,
    labelList& hullEdges,
    labelList& hullFaces,
    labelList& hullCells,
    labelListList& ringEntities
)
{
    // [1] hullEdges is an ordered list of edge-labels around eIndex,
    //     but not connected to it.
    //      - Ordering is in the same manner as vertexHull.
    // [2] hullFaces is an ordered list of face-labels connected to eIndex.
    //      - Ordering is in the same manner as vertexHull.
    // [3] hullCells is an ordered list of cell-labels connected to eIndex.
    //      - For boundary hulls, the last cell label is -1
    // [4] ringEntities are edges and faces connected to eIndex[0] and eIndex[1]
    //      - ringEntities[0]: edges connected to eIndex[0]
    //      - ringEntities[1]: faces connected to eIndex[0]
    //      - ringEntities[2]: edges connected to eIndex[1]
    //      - ringEntities[3]: faces connected to eIndex[1]

    bool found;
    label otherPoint = -1, nextPoint = -1;

    // Obtain a reference to this edge, and its edgeFaces
    const edge& edgeToCheck = edges[eIndex];
    const labelList& eFaces = edgeFaces[eIndex];

    // Loop through all faces of this edge and add them to hullFaces
    forAll(eFaces, faceI)
    {
        const face& faceToCheck = faces[eFaces[faceI]];

        // Find the isolated point on this face,
        // and compare it with vertexHull
        meshOps::findIsolatedPoint
        (
            faceToCheck,
            edgeToCheck,
            otherPoint,
            nextPoint
        );

        found = false;

        forAll(vertexHull, indexI)
        {
            if (vertexHull[indexI] == otherPoint)
            {
                // Fill in the position of this face on the hull
                hullFaces[indexI] = eFaces[faceI];

                // Obtain edges connected to top and bottom
                // vertices of edgeToCheck
                const labelList& fEdges = faceEdges[hullFaces[indexI]];

                forAll(fEdges, edgeI)
                {
                    if
                    (
                        edges[fEdges[edgeI]]
                     == edge(edgeToCheck[0], otherPoint)
                    )
                    {
                        ringEntities[0][indexI] = fEdges[edgeI];
                    }

                    if
                    (
                        edges[fEdges[edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        ringEntities[2][indexI] = fEdges[edgeI];
                    }
                }

                // Depending on the orientation of this face,
                // fill in hull cell indices as well
                if (nextPoint == edgeToCheck[0])
                {
                    hullCells[indexI] = owner[eFaces[faceI]];
                }
                else
                if (nextPoint == edgeToCheck[1])
                {
                    hullCells[indexI] = neighbour[eFaces[faceI]];
                }
                else
                {
                    // Something's terribly wrong
                    FatalErrorIn("void meshOps::constructHull()")
                        << nl << " Failed to construct hull. "
                        << nl << " Possibly not a tetrahedral mesh. "
                        << abort(FatalError);
                }

                if (hullCells[indexI] != -1)
                {
                    label nextI = vertexHull.fcIndex(indexI);
                    label nextHullPoint = vertexHull[nextI];
                    const cell& currCell = cells[hullCells[indexI]];

                    // Look for the ring-faces
                    forAll(currCell, faceI)
                    {
                        const face& cFace = faces[currCell[faceI]];

                        // Check if this face contains edgeToCheck[0]
                        if
                        (
                            triFace::compare
                            (
                                triFace(cFace),
                                triFace
                                (
                                    edgeToCheck[0],
                                    otherPoint,
                                    nextHullPoint
                                )
                            )
                        )
                        {
                            ringEntities[1][indexI] = currCell[faceI];
                        }

                        // Check if this face contains edgeToCheck[1]
                        if
                        (
                            triFace::compare
                            (
                                triFace(cFace),
                                triFace
                                (
                                    edgeToCheck[1],
                                    nextHullPoint,
                                    otherPoint
                                )
                            )
                        )
                        {
                            ringEntities[3][indexI] = currCell[faceI];
                        }
                    }

                    // Scan one the faces for the ring-edge
                    const labelList& rFaceEdges =
                    (
                        faceEdges[ringEntities[1][indexI]]
                    );

                    forAll(rFaceEdges, edgeI)
                    {
                        if
                        (
                            edges[rFaceEdges[edgeI]]
                         == edge(otherPoint,nextHullPoint)
                        )
                        {
                            hullEdges[indexI] = rFaceEdges[edgeI];
                            break;
                        }
                    }
                }

                // Done with this index. Break out.
                found = true;
                break;
            }
        }

        // Throw an error if the point wasn't found
        if (!found)
        {
            // Something's terribly wrong
            FatalErrorIn("void meshOps::constructHull()")
                << " Failed to construct hull. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << " Edge: " << eIndex << ":: " << edgeToCheck << nl
                << " edgeFaces: " << eFaces << nl
                << " vertexHull: " << vertexHull
                << abort(FatalError);
        }
    }
}


// Given a set of points and edges, find the shortest path
// between the start and end point, using Dijkstra's algorithm.
//  - Takes a Map of points and edges that use those points.
//  - Edge weights are currently edge-lengths, but can easily be adapted.
//  - Returns true if the endPoint was found by the algorithm.
//  - The Map 'pi' returns a preceding point for every point in 'points'.
//
//  Algorithm is inspired by:
//    Renaud Waldura
//    Dijkstra's Shortest Path Algorithm in Java
//    http://renaud.waldura.com/
inline bool Dijkstra
(
    const Map<point>& points,
    const Map<edge>& edges,
    const label startPoint,
    const label endPoint,
    Map<label>& pi
)
{
    bool foundEndPoint = false;

    // Set of unvisited (Q) / visited (S) points and distances (d)
    labelHashSet Q, S;
    Map<scalar> d;

    // Initialize distances to large values
    forAllConstIter(Map<point>, points, pIter)
    {
        d.insert(pIter.key(), GREAT);
    }

    // Invert edges to make a local pointEdges list
    Map<labelList> localPointEdges;

    forAllConstIter(Map<edge>, edges, eIter)
    {
        const edge& edgeToCheck = eIter();

        forAll(edgeToCheck, pointI)
        {
            if (!localPointEdges.found(edgeToCheck[pointI]))
            {
                localPointEdges.insert(edgeToCheck[pointI], labelList(0));
            }

            meshOps::sizeUpList
            (
                eIter.key(),
                localPointEdges[edgeToCheck[pointI]]
            );
        }
    }

    // Mark the startPoint as having the smallest distance
    d[startPoint] = 0.0;

    // Add the startPoint to the list of unvisited points
    Q.insert(startPoint);

    while (Q.size())
    {
        // Step 1: Find the node with the smallest distance from the start.
        labelHashSet::iterator smallest = Q.begin();

        for
        (
            labelHashSet::iterator iter = ++Q.begin();
            iter != Q.end();
            iter++
        )
        {
            if (d[iter.key()] < d[smallest.key()])
            {
                smallest = iter;
            }
        }

        label pointIndex = smallest.key();
        scalar smallestDistance = d[pointIndex];

        // Move to the visited points list
        S.insert(pointIndex);
        Q.erase(pointIndex);

        // Step 2: Build a list of points adjacent to pointIndex
        //         but not in the visited list
        dynamicLabelList adjacentPoints(10);

        const labelList& pEdges = localPointEdges[pointIndex];

        forAll(pEdges, edgeI)
        {
            const edge& edgeToCheck = edges[pEdges[edgeI]];

            label otherPoint = edgeToCheck.otherVertex(pointIndex);

            if (!S.found(otherPoint))
            {
                adjacentPoints.append(otherPoint);
            }
        }

        // Step 3: Perform distance-based checks for adjacent points
        forAll(adjacentPoints, pointI)
        {
            label adjPoint = adjacentPoints[pointI];

            scalar distance =
            (
                mag(points[adjPoint] - points[pointIndex])
              + smallestDistance
            );

            // Check if the end-point has been touched.
            if (adjPoint == endPoint)
            {
                foundEndPoint = true;
            }

            if (distance < d[adjPoint])
            {
                // Update to the shorter distance
                d[adjPoint] = distance;

                // Update the predecessor
                if (pi.found(adjPoint))
                {
                    pi[adjPoint] = pointIndex;
                }
                else
                {
                    pi.insert(adjPoint, pointIndex);
                }

                // Add to the list of unvisited points
                Q.insert(adjPoint);
            }
        }
    }

    /*
    // Write out the path
    if (debug > 3)
    {
        if (foundEndPoint)
        {
            dynamicLabelList pathNodes(50);

            label currentPoint = endPoint;

            while (currentPoint != startPoint)
            {
                pathNodes.append(currentPoint);

                currentPoint = pi[currentPoint];
            }

            pathNodes.append(startPoint);

            pathNodes.shrink();

            writeVTK
            (
                "DijkstraPath_"
              + Foam::name(startPoint)
              + '_'
              + Foam::name(endPoint),
                pathNodes,
                0
            );
        }
    }
    */

    return foundEndPoint;
}


// Select a list of elements from connectivity,
// and output to a VTK format
inline void writeVTK
(
    const polyMesh& mesh,
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const UList<point>& meshPoints,
    const UList<edge>& edges,
    const UList<face>& faces,
    const UList<cell>& cells,
    const UList<label>& owner,
    const UList<scalar>& scalField,
    const UList<label>& lablField,
    const UList<vector>& vectField
)
{
    label nTotalCells = 0;
    label nPoints = 0, nCells = 0;

    // Estimate a size for points and cellPoints
    pointField points(6*cList.size());

    // Connectivity lists
    labelListList cpList(cList.size());

    // Create a map for local points
    Map<label> pointMap, reversePointMap, reverseCellMap;

    forAll(cList, cellI)
    {
        if (cList[cellI] < 0)
        {
            continue;
        }

        bool isPolyhedron = false;

        switch (primitiveType)
        {
            // Are we looking at points?
            case 0:
            {
                // Size the list
                cpList[nCells].setSize(1);

                cpList[nCells] = cList[cellI];

                nTotalCells++;

                break;
            }

            // Are we looking at edges?
            case 1:
            {
                // Size the list
                cpList[nCells].setSize(2);

                const edge& tEdge = edges[cList[cellI]];

                cpList[nCells][0] = tEdge[0];
                cpList[nCells][1] = tEdge[1];

                nTotalCells += 2;

                break;
            }

            // Are we looking at faces?
            case 2:
            {
                const face& tFace = faces[cList[cellI]];

                // Size the list and assign
                cpList[nCells] = tFace;

                nTotalCells += tFace.size();

                break;
            }

            // Are we looking at cells?
            case 3:
            {
                const cell& tCell = cells[cList[cellI]];

                // Is face information provided?
                if (faces.size())
                {
                    if (tCell.size() == 4)
                    {
                        // Point-ordering for tetrahedra
                        const face& baseFace = faces[tCell[0]];
                        const face& checkFace = faces[tCell[1]];

                        // Size the list
                        cpList[nCells].setSize(4);

                        // Get the fourth point
                        label apexPoint =
                        (
                            meshOps::findIsolatedPoint(baseFace, checkFace)
                        );

                        // Something's wrong with connectivity.
                        if (apexPoint == -1)
                        {
                            FatalErrorIn
                            (
                                "void meshOps::writeVTK\n"
                                "(\n"
                                "    const polyMesh& mesh,\n"
                                "    const word& name,\n"
                                "    const labelList& cList,\n"
                                "    const label primitiveType,\n"
                                "    const UList<point>& meshPoints,\n"
                                "    const UList<edge>& edges,\n"
                                "    const UList<face>& faces,\n"
                                "    const UList<cell>& cells,\n"
                                "    const UList<label>& owner\n"
                                "    const UList<scalar>& scalField,\n"
                                "    const UList<label>& lablField,\n"
                                "    const UList<vector>& vectField\n"
                                ")\n"
                            )
                                << "Cell: " << cList[cellI]
                                << ":: " << tCell
                                << " has inconsistent connectivity."
                                << abort(FatalError);
                        }

                        // Write-out in order
                        label ownCell = owner[tCell[0]];

                        if (ownCell == cList[cellI])
                        {
                            cpList[nCells][0] = baseFace[2];
                            cpList[nCells][1] = baseFace[1];
                            cpList[nCells][2] = baseFace[0];
                            cpList[nCells][3] = apexPoint;
                        }
                        else
                        {
                            cpList[nCells][0] = baseFace[0];
                            cpList[nCells][1] = baseFace[1];
                            cpList[nCells][2] = baseFace[2];
                            cpList[nCells][3] = apexPoint;
                        }

                        nTotalCells += 4;
                    }
                    else
                    if (tCell.size() > 4)
                    {
                        // Point ordering for polyhedra
                        isPolyhedron = true;

                        // First obtain the face count
                        label npF = 0;

                        forAll(tCell, faceI)
                        {
                            npF += faces[tCell[faceI]].size();
                        }

                        // Size the list
                        cpList[nCells].setSize(tCell.size() + npF + 1);

                        // Fill in facePoints
                        label nP = 0;

                        // Fill in the number of faces
                        cpList[nCells][nP++] = tCell.size();

                        forAll(tCell, faceI)
                        {
                            const face& checkFace = faces[tCell[faceI]];

                            if (checkFace.empty())
                            {
                                FatalErrorIn("void meshOps::writeVTK()")
                                    << " Empty face: " << tCell[faceI]
                                    << abort(FatalError);
                            }

                            // First fill in face size
                            cpList[nCells][nP++] = checkFace.size();

                            // Next fill in points in order
                            if (owner[tCell[faceI]] == cList[cellI])
                            {
                                forAll(checkFace, pI)
                                {
                                    cpList[nCells][nP++] = checkFace[pI];
                                }
                            }
                            else
                            {
                                forAllReverse(checkFace, pI)
                                {
                                    cpList[nCells][nP++] = checkFace[pI];
                                }
                            }
                        }

                        nTotalCells += (tCell.size() + npF + 1);
                    }
                    else
                    {
                        FatalErrorIn("void meshOps::writeVTK()")
                            << " Wrong cell size: " << tCell << nl
                            << " Index: " << cList[cellI] << nl
                            << abort(FatalError);
                    }
                }
                else
                {
                    // No face information.
                    // Assume cell-to-node information.
                    if (tCell.size() == 4)
                    {
                        // Build a face out of first three points.
                        triPointRef tpr
                        (
                            meshPoints[tCell[0]],
                            meshPoints[tCell[1]],
                            meshPoints[tCell[2]]
                        );

                        // Fetch fourth point
                        const point& d = meshPoints[tCell[3]];

                        scalar dDotN = ((d - tpr.a()) & tpr.normal());

                        // Size the list
                        cpList[nCells].setSize(4);

                        if (dDotN > 0.0)
                        {
                            // Correct orientation
                            cpList[nCells][0] = tCell[0];
                            cpList[nCells][1] = tCell[1];
                            cpList[nCells][2] = tCell[2];
                            cpList[nCells][3] = tCell[3];
                        }
                        else
                        {
                            // Flip triangle
                            cpList[nCells][0] = tCell[2];
                            cpList[nCells][1] = tCell[1];
                            cpList[nCells][2] = tCell[0];
                            cpList[nCells][3] = tCell[3];
                        }

                        nTotalCells += 4;
                    }
                    else
                    {
                        FatalErrorIn("void meshOps::writeVTK()")
                            << " Wrong cell size: " << tCell << nl
                            << " Index: " << cList[cellI] << nl
                            << abort(FatalError);
                    }
                }

                break;
            }

            default:
            {
                FatalErrorIn("void meshOps::writeVTK() const")
                    << " Incorrect primitiveType: "
                    << primitiveType
                    << abort(FatalError);
            }
        }

        // Renumber to local ordering
        if (isPolyhedron)
        {
            // Polyhedra also contain face sizes, so use a different scheme
            label pCtr = 0;

            // Fetch number of faces
            label npF = cpList[nCells][pCtr++];

            for (label i = 0; i < npF; i++)
            {
                // Fetch number of points
                label npP = cpList[nCells][pCtr++];

                // Read points in order
                for (label j = 0; j < npP; j++)
                {
                    // Check if this point was added to the map
                    if (!pointMap.found(cpList[nCells][pCtr]))
                    {
                        // Resize pointField if necessary
                        if (nPoints >= points.size())
                        {
                            points.setSize(2 * nPoints);
                        }

                        // Point was not found, so add it
                        points[nPoints] = meshPoints[cpList[nCells][pCtr]];

                        // Update the map
                        pointMap.insert(cpList[nCells][pCtr], nPoints);
                        reversePointMap.insert(nPoints, cpList[nCells][pCtr]);

                        // Increment the number of points
                        nPoints++;
                    }

                    // Renumber it.
                    cpList[nCells][pCtr] = pointMap[cpList[nCells][pCtr]];

                    // Increment point counter
                    pCtr++;
                }
            }
        }
        else
        {
            forAll(cpList[nCells], pointI)
            {
                // Check if this point was added to the map
                if (!pointMap.found(cpList[nCells][pointI]))
                {
                    // Resize pointField if necessary
                    if (nPoints >= points.size())
                    {
                        points.setSize(2 * nPoints);
                    }

                    // Point was not found, so add it
                    points[nPoints] = meshPoints[cpList[nCells][pointI]];

                    // Update the map
                    pointMap.insert(cpList[nCells][pointI], nPoints);
                    reversePointMap.insert(nPoints, cpList[nCells][pointI]);

                    // Increment the number of points
                    nPoints++;
                }

                // Renumber it.
                cpList[nCells][pointI] = pointMap[cpList[nCells][pointI]];
            }
        }

        // Update the cell map.
        reverseCellMap.insert(nCells, cList[cellI]);

        nCells++;
    }

    // Finally write it out
    meshOps::writeVTK
    (
        mesh,
        name,
        nPoints,
        nCells,
        nTotalCells,
        points,
        cpList,
        primitiveType,
        reversePointMap,
        reverseCellMap,
        scalField,
        lablField,
        vectField
    );
}


// Actual routine to write out the VTK file
inline void writeVTK
(
    const polyMesh& mesh,
    const word& name,
    const label nPoints,
    const label nCells,
    const label nTotalCells,
    const vectorField& points,
    const labelListList& cpList,
    const label primitiveType,
    const Map<label>& reversePointMap,
    const Map<label>& reverseCellMap,
    const UList<scalar>& scalField,
    const UList<label>& lablField,
    const UList<vector>& vectField
)
{
    // Make the directory
    fileName dirName(mesh.time().path()/"VTK"/mesh.time().timeName());

    mkDir(dirName);

    // Open stream for output
    OFstream file(dirName/name+".vtk");

    // Write out the header
    file << "# vtk DataFile Version 2.0" << nl
         << name << ".vtk" << nl
         << "ASCII" << nl
         << "DATASET UNSTRUCTURED_GRID" << nl
         << "POINTS " << nPoints << " double" << nl;

    for (label i = 0; i < nPoints; i++)
    {
        file << setprecision(15)
             << points[i].x() << ' '
             << points[i].y() << ' '
             << points[i].z() << ' '
             << nl;
    }

    file << "CELLS " << nCells << " " << nTotalCells + nCells << endl;

    if (cpList.size())
    {
        forAll(cpList, i)
        {
            if (cpList[i].size())
            {
                file << cpList[i].size() << ' ';

                forAll(cpList[i], j)
                {
                    file << cpList[i][j] << ' ';
                }

                file << nl;
            }
        }
    }
    else
    {
        // List of points
        for (label i = 0; i < nPoints; i++)
        {
            file << 1 << ' ' << i << nl;
        }
    }

    file << "CELL_TYPES " << nCells << endl;

    if (cpList.size())
    {
        forAll(cpList, i)
        {
            switch (primitiveType)
            {
                // Vertex
                case 0:
                {
                    file << "1" << nl;
                    break;
                }

                // Edge
                case 1:
                {
                    file << "3" << nl;
                    break;
                }

                // General Polygonal Face
                case 2:
                {
                    file << "7" << nl;
                    break;
                }

                // Cell
                case 3:
                {
                    if (cpList[i].size() == 4)
                    {
                        // Tetrahedron
                        file << "10" << nl;
                    }
                    else
                    if (cpList[i].size() > 4)
                    {
                        // Polyhedron
                        file << "42" << nl;
                    }

                    break;
                }

                default:
                {
                    FatalErrorIn("void meshOps::writeVTK() const")
                        << " Incorrect primitiveType: "
                        << primitiveType
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        // List of points
        for (label i = 0; i < nPoints; i++)
        {
            // Vertex
            file << '1' << nl;
        }
    }

    label nCellFields = 0, nPointFields = 0;

    if (reverseCellMap.size())
    {
        nCellFields++;
    }

    if (scalField.size() == nCells && scalField.size())
    {
        nCellFields++;
    }

    if (lablField.size() == nCells && lablField.size())
    {
        nCellFields++;
    }

    if (vectField.size() == nCells && vectField.size())
    {
        nCellFields++;
    }

    // Write out indices for visualization.
    if (nCellFields)
    {
        file << "CELL_DATA " << nCells << endl;

        file << "FIELD CellFields " << nCellFields << endl;

        if (reverseCellMap.size())
        {
            file << "CellIds 1 " << nCells << " int" << endl;

            for (label i = 0; i < nCells; i++)
            {
                file << reverseCellMap[i] << ' ';
            }

            file << endl;
        }

        if (scalField.size() == nCells)
        {
            file << "CellScalars 1 " << nCells << " double" << endl;

            forAll(scalField, i)
            {
                file << scalField[i] << ' ';
            }

            file << endl;
        }

        if (lablField.size() == nCells)
        {
            file << "CellLabels 1 " << nCells << " int" << endl;

            forAll(lablField, i)
            {
                file << lablField[i] << ' ';
            }

            file << endl;
        }

        if (vectField.size() == nCells)
        {
            file << "CellVectors 3 " << nCells << " double" << endl;

            forAll(vectField, i)
            {
                file << vectField[i].x() << ' '
                     << vectField[i].y() << ' '
                     << vectField[i].z() << ' ';
            }

            file << endl;
        }
    }

    if (reversePointMap.size())
    {
        nPointFields++;
    }

    if (scalField.size() == nPoints && scalField.size())
    {
        nPointFields++;
    }

    if (lablField.size() == nPoints && lablField.size())
    {
        nPointFields++;
    }

    if (vectField.size() == nPoints && vectField.size())
    {
        nPointFields++;
    }

    // Write out indices for visualization.
    if (nPointFields)
    {
        file << "POINT_DATA " << nPoints << endl;

        file << "FIELD PointFields " << nPointFields << endl;

        if (reversePointMap.size())
        {
            file << "PointIds 1 " << nPoints << " int" << endl;

            for (label i = 0; i < nPoints; i++)
            {
                file << reversePointMap[i] << ' ';
            }

            file << endl;
        }

        if (scalField.size() == nPoints)
        {
            file << "PointScalars 1 " << nPoints << " double" << endl;

            forAll(scalField, i)
            {
                file << scalField[i] << ' ';
            }

            file << endl;
        }

        if (lablField.size() == nPoints)
        {
            file << "PointLabels 1 " << nPoints << " int" << endl;

            forAll(lablField, i)
            {
                file << lablField[i] << ' ';
            }

            file << endl;
        }

        if (vectField.size() == nPoints)
        {
            file << "PointVectors 3 " << nPoints << " double" << endl;

            forAll(vectField, i)
            {
                file << vectField[i].x() << ' '
                     << vectField[i].y() << ' '
                     << vectField[i].z() << ' ';
            }

            file << endl;
        }
    }
}


} // End namespace meshOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
