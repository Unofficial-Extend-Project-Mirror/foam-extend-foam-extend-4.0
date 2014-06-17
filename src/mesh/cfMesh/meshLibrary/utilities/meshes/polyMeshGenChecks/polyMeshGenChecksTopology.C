/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenChecks.H"
#include "polyMeshGenAddressing.H"
#include "cell.H"
#include "Map.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace polyMeshGenChecks
{

bool checkPoints
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    label nFaceErrors = 0;
    label nCellErrors = 0;

    const VRWGraph& pf = mesh.addressingData().pointFaces();

    forAll(pf, pointI)
    {
        if( pf.sizeOfRow(pointI) == 0 )
        {
            WarningIn
            (
                "bool checkPoints"
                "(const polyMeshGen&, const bool, labelHashSet*)"
            )   << "Point " << pointI << " not used by any faces." << endl;

            if( setPtr )
                setPtr->insert(pointI);

            ++nFaceErrors;
        }
    }

    const VRWGraph& pc = mesh.addressingData().pointCells();

    forAll(pc, pointI)
    {
        if( pc.sizeOfRow(pointI) == 0 )
        {
            WarningIn
            (
                "bool checkPoints"
                "(const polyMeshGen&, const bool, labelHashSet*)"
            )   << "Point " << pointI << " not used by any cells." << endl;

            if( setPtr )
                setPtr->insert(pointI);

            ++nCellErrors;
        }
    }

    reduce(nFaceErrors, sumOp<label>());
    reduce(nCellErrors, sumOp<label>());

    if( nFaceErrors > 0 || nCellErrors > 0 )
    {
        WarningIn
        (
            "bool checkPoints"
            "(const polyMeshGen&, const bool, labelHashSet*)"
        )   << "Error in point usage detected: " << nFaceErrors
            << " unused points found in the mesh.  This mesh is invalid."
            << endl;

        return true;
    }
    else
    {
        if( report )
            Info << "Point usage check OK.\n" << endl;

        return false;
    }
}

bool checkUpperTriangular
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    // Check whether internal faces are ordered in the upper triangular order
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const cellListPMG& cells = mesh.cells();
    const VRWGraph& cc = mesh.addressingData().cellCells();

    const label internal = mesh.nInternalFaces();

    labelList checkInternalFaces(internal, -1);

    label nChecks = 0;

    bool error = false;

    // Loop through faceCells once more and make sure that for internal cell
    // the first label is smaller
    for(label faceI=0;faceI<internal;++faceI)
    {
        if( own[faceI] >= nei[faceI] )
        {
            if( report )
            {
                Pout<< "bool checkUpperTriangular(const polyMeshGen&, "
                    << "const bool, labelHashSet*) : " << endl
                    << "face " << faceI
                    << " has the owner label greater than neighbour:" << endl
                    << own[faceI] << tab << nei[faceI] << endl;
            }

            if( setPtr )
                setPtr->insert(faceI);

            error  = true;
        }
    }

    // Loop through all cells. For each cell, find the face that is internal and
    // add it to the check list (upper triangular order).
    // Once the list is completed, check it against the faceCell list
    forAll(cells, cellI)
    {
        const labelList& curFaces = cells[cellI];

        // Using the fact that cell neighbour always appear
        // in the increasing order
        boolList usedNbr(cc.sizeOfRow(cellI), false);

        for(label nSweeps=0;nSweeps<usedNbr.size();++nSweeps)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = cells.size();

            forAllRow(cc, cellI, nbrI)
            {
                const label neiI = cc(cellI, nbrI);
                if( (neiI > cellI) && !usedNbr[nbrI] && (neiI < minNei) )
                {
                    nextNei = nbrI;
                    minNei = neiI;
                }
            }

            if( nextNei > -1 )
            {
                // Mark this neighbour as used
                usedNbr[nextNei] = true;

                forAll(curFaces, faceI)
                {
                    if( curFaces[faceI] < internal )
                    {
                        if( nei[curFaces[faceI]] == cc(cellI, nextNei) )
                        {
                            checkInternalFaces[nChecks] = curFaces[faceI];
                            ++nChecks;

                            break;
                        }
                    }
                }
            }
        }
    }

    // Check list created. If everything is OK, the face label is equal to index
    forAll(checkInternalFaces, faceI)
    {
        if( checkInternalFaces[faceI] != faceI )
        {
            error = true;

            Pout<< "bool checkUpperTriangular(const polyMeshGen&, const bool"
                << ", labelHashSet*) : " << endl
                << "face " << faceI << " out of position. Markup label: "
                << checkInternalFaces[faceI] << ". All subsequent faces will "
                << "also be out of position. Please check the mesh manually."
                << endl;

            if( setPtr )
                setPtr->insert(faceI);

            break;
        }
    }

    reduce(error, orOp<bool>());

    if( error )
    {
        WarningIn
        (
            "bool checkUpperTriangular(const polyMeshGen&, const bool"
            ", labelHashSet*)"
        )   << "Error in face ordering: faces not in upper triangular order!"
            << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Upper triangular ordering OK.\n" << endl;

        return false;
    }
}

bool checkCellsZipUp
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    label nOpenCells = 0;

    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided) reduction(+ : nOpenCells)
    # endif
    forAll(cells, cellI)
    {
        const labelList& c = cells[cellI];

        DynList<edge> cellEdges;
        DynList<label> edgeUsage;

        forAll(c, faceI)
        {
            const face& f = faces[c[faceI]];

            forAll(f, eI)
            {
                const edge e = f.faceEdge(eI);

                const label pos = cellEdges.containsAtPosition(e);

                if( pos < 0 )
                {
                    cellEdges.append(e);
                    edgeUsage.append(1);
                }
                else
                {
                    ++edgeUsage[pos];
                }
            }
        }

        DynList<edge> singleEdges;

        forAll(edgeUsage, edgeI)
        {
            if( edgeUsage[edgeI] == 1 )
            {
                singleEdges.append(cellEdges[edgeI]);
            }
            else if( edgeUsage[edgeI] != 2 )
            {
                WarningIn
                (
                    "bool checkCellsZipUp(const polyMeshGen&,"
                    "const bool, labelHashSet*)"
                )   << "edge " << cellEdges[edgeI] << " in cell " << cellI
                    << " used " << edgeUsage[edgeI] << " times. " << endl
                    << "Should be 1 or 2 - serious error in mesh structure"
                    << endl;

                if( setPtr )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    setPtr->insert(cellI);
                }
            }
        }

        if( singleEdges.size() > 0 )
        {
            if( report )
            {
                Pout<< "bool checkCellsZipUp(const polyMeshGen&, const bool"
                    << ", labelHashSet*) : " << endl
                    << "Cell " << cellI << " has got " << singleEdges.size()
                    << " unmatched edges: " << singleEdges << endl;
            }

            if( setPtr )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                setPtr->insert(cellI);
            }

            ++nOpenCells;
        }
    }

    reduce(nOpenCells, sumOp<label>());

    if( nOpenCells > 0 )
    {
        WarningIn
        (
            "bool checkCellsZipUp(const polyMeshGen&,"
            " const bool, labelHashSet*)"
        )   << nOpenCells
            << " open cells found.  Please use the mesh zip-up tool. "
            << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Topological cell zip-up check OK.\n" << endl;

        return false;
    }
}

// Vertices of face within point range and unique.
bool checkFaceVertices
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    // Check that all vertex labels are valid
    const faceListPMG& faces = mesh.faces();

    label nErrorFaces = 0;
    const label nPoints = mesh.points().size();

    forAll(faces, fI)
    {
        const face& curFace = faces[fI];

        if( min(curFace) < 0 || max(curFace) > nPoints )
        {
            WarningIn
            (
                "bool checkFaceVertices("
                "const polyMesgGen&, const bool, labelHashSet*)"
            )   << "Face " << fI << " contains vertex labels out of range: "
                << curFace << " Max point index = " << nPoints-1 << endl;

            if( setPtr )
                setPtr->insert(fI);

            ++nErrorFaces;
        }

        // Uniqueness of vertices
        labelHashSet facePoints(2*curFace.size());

        forAll(curFace, fp)
        {
            bool inserted = facePoints.insert(curFace[fp]);

            if( !inserted )
            {
                WarningIn
                (
                    "bool checkFaceVertices("
                    "const polyMeshGen&, const bool, labelHashSet*)"
                )   << "Face " << fI << " contains duplicate vertex labels: "
                    << curFace << endl;

                if( setPtr )
                    setPtr->insert(fI);

                ++nErrorFaces;
            }
        }
    }

    reduce(nErrorFaces, sumOp<label>());

    if( nErrorFaces > 0 )
    {
        SeriousErrorIn
        (
            "bool checkFaceVertices("
            "const polyMeshGen&, const bool, labelHashSet*)"
        )   << "const bool, labelHashSet*) const: "
            << nErrorFaces << " faces with invalid vertex labels found"
            << endl;

        return true;
    }
    else
    {
        if( report )
            Info<< "Face vertices OK.\n" << endl;

        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace polyMeshGenChecks

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
