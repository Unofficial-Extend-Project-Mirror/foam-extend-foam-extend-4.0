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

Description
    Refines a hex mesh by 2x2x2 cell splitting.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "pointMesh.H"
#include "argList.H"
#include "foamTime.H"
#include "hexRef8.H"
#include "cellSet.H"
#include "OFstream.H"
#include "meshTools.H"
#include "IFstream.H"
#include "directTopoChange.H"
#include "mapPolyMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");
    argList::validArgs.append("cellSet");
#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const pointMesh& pMesh = pointMesh::New(mesh);

    word cellSetName(args.args()[1]);
    bool overwrite = args.optionFound("overwrite");

    Info<< "Reading cells to refine from cellSet " << cellSetName
        << nl << endl;

    cellSet cellsToRefine(mesh, cellSetName);

    Info<< "Read " << returnReduce(cellsToRefine.size(), sumOp<label>())
        << " cells to refine from cellSet " << cellSetName << nl
        << endl;


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);

    // Read point fields
    PtrList<pointScalarField> psFlds;
    ReadFields(pMesh, objects, psFlds);

    PtrList<pointVectorField> pvFlds;
    ReadFields(pMesh, objects, pvFlds);



    // Construct refiner without unrefinement. Read existing point/cell level.
    hexRef8 meshCutter(mesh);

    // Some stats
    Info<< "Read mesh:" << nl
        << "    cells:" << mesh.globalData().nTotalCells() << nl
        << "    faces:" << mesh.globalData().nTotalFaces() << nl
        << "    points:" << mesh.globalData().nTotalPoints() << nl
        << "    cellLevel :"
        << " min:" << gMin(meshCutter.cellLevel())
        << " max:" << gMax(meshCutter.cellLevel()) << nl
        << "    pointLevel :"
        << " min:" << gMin(meshCutter.pointLevel())
        << " max:" << gMax(meshCutter.pointLevel()) << nl
        << endl;


    // Maintain 2:1 ratio
    labelList newCellsToRefine
    (
        meshCutter.consistentRefinement
        (
            cellsToRefine.toc(),
            true                  // extend set
        )
    );

    // Mesh changing engine.
    directTopoChange meshMod(mesh);

    // Play refinement commands into mesh changer.
    meshCutter.setRefinement(newCellsToRefine, meshMod);

    if (!overwrite)
    {
        runTime++;
    }

    // Create mesh, return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Update numbering of cells/vertices.
    meshCutter.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    Pout<< "Refined from " << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh.globalData().nTotalCells() << " cells."
        << nl << endl;

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
        meshCutter.setInstance(oldInstance);
    }
    Info<< "Writing mesh to " << runTime.timeName() << endl;

    mesh.write();
    meshCutter.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
