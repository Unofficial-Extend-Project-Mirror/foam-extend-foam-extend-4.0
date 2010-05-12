#include "printMeshStats.H"
#include "polyMesh.H"
#include "globalMeshData.H"

#include "hexMatcher.H"
#include "wedgeMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetWedgeMatcher.H"
#include "tetMatcher.H"


void Foam::printMeshStats(const polyMesh& mesh, const bool allTopology)
{
    Pout<< "Mesh stats" << nl
        << "    points:           " << mesh.points().size() << nl;

    if (mesh.nInternalPoints() != -1)
    {
        Pout<< "    internal points:  " << mesh.nInternalPoints() << nl;
    }

    if (allTopology && mesh.nInternalPoints() != -1)
    {
        Pout<< "    edges:            " << mesh.nEdges() << nl
            << "    internal edges:   " << mesh.nInternalEdges() << nl
            << "    internal edges using one boundary point:   "
            << mesh.nInternal1Edges()-mesh.nInternal0Edges() << nl
            << "    internal edges using two boundary points:  "
            << mesh.nInternalEdges()-mesh.nInternal1Edges() << nl;
    }

    Pout<< "    faces:            " << mesh.faces().size() << nl
        << "    internal faces:   " << mesh.faceNeighbour().size() << nl
        << "    cells:            " << mesh.cells().size() << nl
        << "    boundary patches: " << mesh.boundaryMesh().size() << nl
        << "    point zones:      " << mesh.pointZones().size() << nl
        << "    face zones:       " << mesh.faceZones().size() << nl
        << "    cell zones:       " << mesh.cellZones().size() << nl
        << endl;

    if (Pstream::parRun())
    {
        const globalMeshData& parData = mesh.globalData();

        Info<< "Overall stats" << nl
            << "    points:   " << parData.nTotalPoints() << nl
            << "    faces:    " << parData.nTotalFaces() << nl
            << "    cells:    " << parData.nTotalCells() << nl
            << endl;
    }

    // Construct shape recognizers
    hexMatcher hex;
    prismMatcher prism;
    wedgeMatcher wedge;
    pyrMatcher pyr;
    tetWedgeMatcher tetWedge;
    tetMatcher tet;

    // Counters for different cell types
    label nHex = 0;
    label nWedge = 0;
    label nPrism = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTetWedge = 0;
    label nUnknown = 0;

    for(label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (hex.isA(mesh, cellI))
        {
            nHex++;
        }
        else if (tet.isA(mesh, cellI))
        {
            nTet++;
        }
        else if (pyr.isA(mesh, cellI))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, cellI))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, cellI))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, cellI))
        {
            nTetWedge++;
        }
        else
        {
            nUnknown++;
        }
    }

    Pout<< "Number of cells of each type: " << nl
        << "    hexahedra:     " << nHex << nl
        << "    prisms:        " << nPrism << nl
        << "    wedges:        " << nWedge << nl
        << "    pyramids:      " << nPyr << nl
        << "    tet wedges:    " << nTetWedge << nl
        << "    tetrahedra:    " << nTet << nl
        << "    polyhedra:     " << nUnknown
        << nl << endl;
}
