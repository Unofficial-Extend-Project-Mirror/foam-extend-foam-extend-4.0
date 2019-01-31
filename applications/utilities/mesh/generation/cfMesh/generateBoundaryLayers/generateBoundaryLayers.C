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
    Generates boundary layers in the existing mesh, based on the settings
    given in meshDict. It also performs necessary quality optimisation.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGenModifier.H"
#include "meshOptimizer.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void generateLayer
(
    polyMeshGen& mesh,
    const dictionary& meshDict,
    const bool layers2D
)
{
    boundaryLayers bl(mesh);

    if( layers2D )
        bl.activate2DMode();

    if( meshDict.found("boundaryLayers") )
    {
        const dictionary& bndLayers = meshDict.subDict("boundaryLayers");

        if( bndLayers.found("nLayers") )
        {
            const label nLayers = readLabel(bndLayers.lookup("nLayers"));

            if( nLayers > 0 )
                bl.addLayerForAllPatches();
        }
        else if( bndLayers.found("patchBoundaryLayers") )
        {
            const dictionary& patchLayers =
                bndLayers.subDict("patchBoundaryLayers");
            const wordList createLayers = patchLayers.toc();

            forAll(createLayers, patchI)
                bl.addLayerForPatch(createLayers[patchI]);
        }
    }
    else
    {
        bl.addLayerForAllPatches();
    }
}

void meshOptimisation(polyMeshGen& mesh)
{
    meshOptimizer mOpt(mesh);

    mOpt.optimizeMeshFV();
    mOpt.optimizeLowQualityFaces();
    mOpt.untangleMeshFV();
    mOpt.optimizeBoundaryLayer();
    mOpt.untangleMeshFV();
}

void layerRefinement(polyMeshGen& mesh, const dictionary& meshDict)
{
    if( meshDict.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh);

        refineBoundaryLayers::readSettings(meshDict, refLayers);

        refLayers.refineLayers();

        meshOptimizer(mesh).untangleBoundaryLayer();
    }
}

int main(int argc, char *argv[])
{
    argList::validOptions.insert("2DLayers", "bool");

#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //- load the mesh from disk
    polyMeshGen pmg(runTime);
    pmg.read();

    bool is2DLayer(false);
    if( args.options().found("2DLayers") )
        is2DLayer = true;

    //- generate the initial boundary layer
    generateLayer(pmg, meshDict, is2DLayer);

    //- optimisation of mesh quality
    meshOptimisation(pmg);

    //- perform layer refinement
    layerRefinement(pmg, meshDict);

    Info << "Writing mesh" << endl;
    pmg.write();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
