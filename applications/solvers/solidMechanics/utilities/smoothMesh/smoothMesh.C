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

Application
    smoothMesh

Description
    Smoothing mesh using Laplacian smoothing

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "boolList.H"
#include "emptyPolyPatch.H"
#include "twoDPointCorrector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    Foam::argList::validOptions.insert("smoothPatches", "");

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    // bool smoothPatches = args.optionFound("smoothPatches");
    // if(smoothPatches)
    //   Info << "Smoothing patches" << nl << endl;

    // Smooth internal points

    const vectorField& oldPoints = mesh.points();

    vectorField newPoints = oldPoints;

    if (mesh.nGeometricD() == 3)
    {
        Info << "3-D mesh" << endl;

        const labelListList& pointEdges = mesh.pointEdges();
        const edgeList& edges = mesh.edges();

        boolList fixedPoints(newPoints.size(), false);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pointI)
            {
                fixedPoints[meshPoints[pointI]] = true;
            }
        }

        scalarField residual(newPoints.size(), 0);
        label counter = 0;
        do
        {
            counter++;

            forAll(newPoints, pointI)
            {
                if (!fixedPoints[pointI])
                {
                    vector curNewPoint = vector::zero;

                    scalar sumW = 0;

                    forAll(pointEdges[pointI], eI)
                    {
                        label curEdgeIndex = pointEdges[pointI][eI];

                        const edge& curEdge = edges[curEdgeIndex];

                        vector d =
                            newPoints[curEdge.otherVertex(pointI)]
                          - newPoints[pointI];

                        scalar w = 1.0;

                        curNewPoint += w*d;

                        sumW += w;
                    }

                    curNewPoint /= sumW;

                    curNewPoint += newPoints[pointI];

                    residual[pointI] = mag(curNewPoint - newPoints[pointI]);

                    newPoints[pointI] = curNewPoint;
                }
            }

        // smoothed patch points independently
        // but we do not smooth points on the boundary of the patch
        // we reset these points
        // problem: we need to smooth feature edges separately
        // if(smoothPatches)
        //   {
        //      forAll(mesh.boundaryMesh(), patchI)
        //        {
        //          const labelList& meshPoints =
        //            mesh.boundaryMesh()[patchI].meshPoints();

        //          forAll(meshPoints, pointI)
        //            {
        //          label pointID = meshPoints[pointI];
        //          vector curNewPoint = vector::zero;
        //          scalar sumW = 0;

        //          forAll(pointEdges[pointID], eI)
        //            {
        //              label curEdgeIndex = pointEdges[pointID][eI];
        //              const edge& curEdge = edges[curEdgeIndex];

        //              // use only boundary points
        //              label otherPointID = curEdge.otherVertex(pointID);
        //              if(fixedPoints[otherPointID])
        //                {
        //              vector d =
        //                newPoints[otherPointID]
        //                - newPoints[pointID];

        //              scalar w = 1.0;
        //              curNewPoint += w*d;
        //              sumW += w;
        //                }
        //            }

        //          curNewPoint /= sumW;
        //          curNewPoint += newPoints[pointID];
        //          residual[pointID] = mag(curNewPoint - newPoints[pointID]);
        //          newPoints[pointID] = curNewPoint;
        //            }

        //          // reset boundary points
        //          const labelList& boundaryPoints =
            // mesh.boundaryMesh()[patchI].boundaryPoints();
        //          forAll(boundaryPoints, bpi)
        //            {
        //          label boundaryPointID = meshPoints[boundaryPoints[bpi]];
        //          newPoints[boundaryPointID] = oldPoints[boundaryPointID];
        //          residual[boundaryPointID] = 0.0;
        //            }
        //        }
        //   }

            residual /= max(mag(newPoints - oldPoints) + SMALL);
        }
        while(max(residual) > 1e-3);

        Info << "Internal points, max residual: " << max(residual)
            << ", num of iterations: " << counter << endl;

        twoDPointCorrector twoDCorrector(mesh);
        twoDCorrector.correctPoints(newPoints);
    }
    else // 2-D
    {
        Info << "2-D mesh" << endl;

        forAll(mesh.boundaryMesh(), patchI)
        {
            if
            (
                mesh.boundaryMesh()[patchI].type()
             == emptyPolyPatch::typeName
            )
            {
                const labelList& bPoints =
                    mesh.boundaryMesh()[patchI].boundaryPoints();

                const vectorField& oldPatchPoints =
                    mesh.boundaryMesh()[patchI].localPoints();

                vectorField newPatchPoints = oldPatchPoints;

                const labelListList& patchPointEdges =
                    mesh.boundaryMesh()[patchI].pointEdges();

                const edgeList& patchEdges =
                    mesh.boundaryMesh()[patchI].edges();

                boolList fixedPoints(newPatchPoints.size(), false);

                forAll(bPoints, pointI)
                {
                    fixedPoints[bPoints[pointI]] = true;
                }

                scalarField residual(newPatchPoints.size(), 0);
                label counter = 0;
                do
                {
                    counter++;

                    forAll(newPatchPoints, pointI)
                    {
                        if (!fixedPoints[pointI])
                        {
                            vector curNewPatchPoint = vector::zero;

                            scalar sumW = 0;

                            forAll(patchPointEdges[pointI], eI)
                            {
                                label curEdgeIndex =
                                    patchPointEdges[pointI][eI];

                                const edge& curEdge = patchEdges[curEdgeIndex];

                                vector d =
                                    newPatchPoints[curEdge.otherVertex(pointI)]
                                  - newPatchPoints[pointI];

                                scalar w = 1.0;

                                curNewPatchPoint += w*d;

                                sumW += w;
                            }

                            curNewPatchPoint /= sumW;

                            curNewPatchPoint += newPatchPoints[pointI];

                            residual[pointI] =
                                mag(curNewPatchPoint - newPatchPoints[pointI]);

                            newPatchPoints[pointI] = curNewPatchPoint;
                        }
                    }

                    residual /=
                        max(mag(newPatchPoints - oldPatchPoints) + SMALL);
                }
                while(max(residual) > 1e-6);

                Info << "Empty points, max residual: " << max(residual)
                    << ", num of iterations: " << counter << endl;

                const labelList& meshPoints =
                    mesh.boundaryMesh()[patchI].meshPoints();

                forAll(meshPoints, pointI)
                {
                    newPoints[meshPoints[pointI]] = newPatchPoints[pointI];
                }
            }
        }
    }


    twoDPointCorrector twoDCorrector(mesh);
    twoDCorrector.correctPoints(newPoints);

    mesh.movePoints(newPoints);

    runTime++;

    runTime.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
