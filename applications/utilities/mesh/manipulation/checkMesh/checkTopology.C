#include "checkTopology.H"
#include "polyMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOmanip.H"

Foam::label Foam::checkTopology
(
    const polyMesh& mesh,
    const bool allTopology,
    const bool allGeometry
)
{
    label noFailedChecks = 0;

    Pout<< "Checking topology..." << endl;

    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);

    // Check if the boundary processor patches are correct
    mesh.boundaryMesh().checkParallelSync(true);

    {
        pointSet points(mesh, "unusedPoints", mesh.nPoints()/100);
        if (mesh.checkPoints(true, &points))
        {
            noFailedChecks++;

            Pout<< "  <<Writing " << points.size()
                << " unused points to set " << points.name() << endl;
            points.write();
        }
    }

    {
        faceSet faces(mesh, "upperTriangularFace", mesh.nFaces()/100);
        if (mesh.checkUpperTriangular(true, &faces))
        {
            noFailedChecks++;
        }
        if (faces.size() > 0)
        {
            Pout<< "  <<Writing " << faces.size()
                << " unordered faces to set " << faces.name() << endl;
            faces.write();
        }
    }

    if (allTopology)
    {
        cellSet cells(mesh, "zipUpCells", mesh.nCells()/100);
        if (mesh.checkCellsZipUp(true, &cells))
        {
            noFailedChecks++;

            Pout<< "  <<Writing " << cells.size()
                << " cells with over used edges to set " << cells.name()
                << endl;
            cells.write();
        }
    }

    {
        faceSet faces(mesh, "outOfRangeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceVertices(true, &faces))
        {
            noFailedChecks++;

            Pout<< "  <<Writing " << faces.size()
                << " faces with out-of-range or duplicate vertices to set "
                << faces.name() << endl;
            faces.write();
        }
    }

    if (allTopology)
    {
        faceSet faces(mesh, "edgeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceFaces(true, &faces))
        {
            noFailedChecks++;

            Pout<< "  <<Writing " << faces.size()
                << " faces with incorrect edges to set " << faces.name()
                << endl;
            faces.write();
        }
    }

    {
        regionSplit rs(mesh);

        if (rs.nRegions() == 1)
        {
            Info<< "    Number of regions: " << rs.nRegions() << " (OK)."
                << endl;
        
        }
        else
        {
            Info<< "   *Number of regions: " << rs.nRegions() << endl;

            Info<< "    The mesh has multiple regions which are not connected "
                   "by any face." << endl
                << "  <<Writing region information to "
                << mesh.time().timeName()/"cellToRegion"
                << endl;

            labelIOList ctr
            (
                IOobject
                (
                    "cellToRegion",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rs
            );
            ctr.write();
        }
    }

    if (!Pstream::parRun())
    {
        Pout<< "\nChecking patch topology for multiply connected surfaces ..."
            << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Non-manifold points
        pointSet points
        (
            mesh,
            "nonManifoldPoints",
            mesh.nPoints()/100
        );

        Pout.setf(ios_base::left);

        Pout<< "    "
            << setw(20) << "Patch"
            << setw(9) << "Faces"
            << setw(9) << "Points"
            << setw(34) << "Surface topology";
        if (allGeometry)
        {
            Pout<< " Bounding box";
        }
        Pout<< endl;

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

                Pout<< "    "
                    << setw(20) << pp.name()
                    << setw(9) << pp.size()
                    << setw(9) << pp.nPoints();


            primitivePatch::surfaceTopo pTyp = pp.surfaceType();

            if (pp.size() == 0)
            {
                Pout<< setw(34) << "ok (empty)";
            }
            else if (pTyp == primitivePatch::MANIFOLD)
            {
                if (pp.checkPointManifold(true, &points))
                {
                    Pout<< setw(34) << "multiply connected (shared point)";
                }
                else
                {
                    Pout<< setw(34) << "ok (closed singly connected)";
                }

                // Add points on non-manifold edges to make set complete
                pp.checkTopology(false, &points);
            }
            else
            {
                pp.checkTopology(false, &points);

                if (pTyp == primitivePatch::OPEN)
                {
                    Pout<< setw(34) << "ok (non-closed singly connected)";
                }
                else
                {
                    Pout<< setw(34) << "multiply connected (shared edge)";
                }
            }

            if (allGeometry)
            {
                const pointField& pts = pp.points();
                const labelList& mp = pp.meshPoints();

                boundBox bb(vector::zero, vector::zero);
                if (returnReduce(mp.size(), sumOp<label>()) > 0)
                {
                    bb.min() = pts[mp[0]];
                    bb.max() = pts[mp[0]];
                    for (label i = 1; i < mp.size(); i++)
                    {
                        bb.min() = min(bb.min(), pts[mp[i]]);
                        bb.max() = max(bb.max(), pts[mp[i]]);
                    }
                    reduce(bb.min(), minOp<vector>());
                    reduce(bb.max(), maxOp<vector>());
                }
                Pout<< ' ' << bb;
            }
            Pout<< endl;
        }

        if (points.size() > 0)
        {
            Pout<< "  <<Writing " << points.size()
                << " conflicting points to set "
                << points.name() << endl;

            points.write();
        }

        //Pout.setf(ios_base::right);
    }

    // Force creation of all addressing if requested.
    // Errors will be reported as required
    if (allTopology)
    {
        mesh.cells();
        mesh.faces();
        mesh.edges();
        mesh.points();
        mesh.faceOwner();
        mesh.faceNeighbour();
        mesh.cellCells();
        mesh.edgeCells();
        mesh.pointCells();
        mesh.edgeFaces();
        mesh.pointFaces();
        mesh.cellEdges();
        mesh.faceEdges();
        mesh.pointEdges();
    }

    return noFailedChecks;
}
