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

\*---------------------------------------------------------------------------*/

#include "leastSquaresVolPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "cyclicPolyPatch.H"
#include "cyclicGgiPolyPatch.H"
#include "cyclicFvPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "Map.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(leastSquaresVolPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void leastSquaresVolPointInterpolation::makePointFaces() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makePointFaces() : "
            << "constructing point boundary faces addressing"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointBndFacesPtr_ || pointProcFacesPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makePointFaces() const"
        )
            << "point boundary faces adressing already exist"
            << abort(FatalError);
    }

    const pointField& points = mesh().points();
    const labelListList& pointFaces = mesh().pointFaces();
    const labelListList& pointPoints = mesh().pointPoints();

    // Allocate storage for addressing
    pointBndFacesPtr_ = new labelListList(points.size());
    labelListList& pointBndFaces = *pointBndFacesPtr_;

    // Allocate storage for addressing
    pointCyclicFacesPtr_ = new labelListList(points.size());
    labelListList& pointCyclicFaces = *pointCyclicFacesPtr_;

    // Allocate storage for addressing
    pointCyclicBndFacesPtr_ = new labelListList(points.size());
    labelListList& pointCyclicBndFaces = *pointCyclicBndFacesPtr_;

    // Allocate storage for addressing
    pointCyclicGgiFacesPtr_ = new labelListList(points.size());
    labelListList& pointCyclicGgiFaces = *pointCyclicGgiFacesPtr_;

    // Allocate storage for addressing
    pointCyclicGgiBndFacesPtr_ = new labelListList(points.size());
    labelListList& pointCyclicGgiBndFaces = *pointCyclicGgiBndFacesPtr_;

    // Allocate storage for addressing
    pointProcFacesPtr_ = new labelListList(points.size());
    labelListList& pointProcFaces = *pointProcFacesPtr_;

    forAll(pointBndFaces, pointI)
    {
        const labelList& curPointFaces = pointFaces[pointI];

        labelHashSet bndFaceSet;
        labelHashSet cyclicFaceSet;
        labelHashSet cyclicGgiFaceSet;
        labelHashSet procFaceSet;

        forAll(curPointFaces, faceI)
        {
            label faceID = curPointFaces[faceI];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            if (patchID != -1)
            {
                if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == cyclicPolyPatch::typeName
                )
                {
                    cyclicFaceSet.insert(faceID);
                }
                else if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == cyclicGgiPolyPatch::typeName
                )
                {
                    cyclicGgiFaceSet.insert(faceID);
                }
                else if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == processorPolyPatch::typeName
                )
                {
                    procFaceSet.insert(faceID);
                }
                else if
                (
                    (
                        mesh().boundaryMesh()[patchID].type()
                     != emptyPolyPatch::typeName
                    )
                 && (
                        mesh().boundaryMesh()[patchID].type()
                     != wedgePolyPatch::typeName
                    )
                )
                {
                    bndFaceSet.insert(faceID);
                }
                else if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == wedgePolyPatch::typeName
                )
                {
                    if (pointAxisEdges().found(pointI))
                    {
                        bndFaceSet.insert(faceID);
                    }
                }
            }
        }

        pointBndFaces[pointI] = bndFaceSet.toc();
        pointCyclicFaces[pointI] = cyclicFaceSet.toc();
        pointCyclicGgiFaces[pointI] = cyclicGgiFaceSet.toc();

        const labelList& glPoints =
            mesh().globalData().sharedPointLabels();

        bool globalPoint(findIndex(glPoints, pointI) != -1);

        if (!globalPoint)
        {
            labelList allPointProcFaces = procFaceSet.toc();

            // Check for duplicate proc faces
            vectorField allCentres(allPointProcFaces.size());

            forAll(allCentres, faceI)
            {
                label faceID = allPointProcFaces[faceI];
                label patchID = mesh().boundaryMesh().whichPatch(faceID);
                label start = mesh().boundaryMesh()[patchID].start();
                label localFaceID = faceID - start;

                allCentres[faceI] =
                    mesh().C().boundaryField()[patchID][localFaceID];
            }

            boundBox bb(vectorField(points, pointPoints[pointI]), false);
            scalar tol = 0.001*mag(bb.max() - bb.min());

            vectorField centres(allCentres.size(), vector::zero);
            pointProcFaces[pointI] = labelList(allCentres.size(), -1);

            label nCentres = 0;
            forAll(allCentres, faceI)
            {
                bool duplicate = false;
                for (label i=0; i<nCentres; i++)
                {
                    if
                    (
                        mag
                        (
                            centres[i]
                          - allCentres[faceI]
                        )
                      < tol
                    )
                    {
                        duplicate = true;
                        break;
                    }
                }

                if (!duplicate)
                {
                    centres[nCentres] = allCentres[faceI];
                    pointProcFaces[pointI][nCentres] =
                        allPointProcFaces[faceI];
                    nCentres++;
                }
            }

            pointProcFaces[pointI].setSize(nCentres);
        }
    }

    // Find cyclic boundary faces
    forAll(pointCyclicBndFaces, pointI)
    {
        if (pointCyclicFaces[pointI].size())
        {
            if (pointBndFaces[pointI].size())
            {
                label faceID = pointCyclicFaces[pointI][0];
                label patchID = mesh().boundaryMesh().whichPatch(faceID);

                label start = mesh().boundaryMesh()[patchID].start();
                label localFaceID = faceID - start;

                const cyclicPolyPatch& cycPatch =
                    refCast<const cyclicPolyPatch>
                    (
                        mesh().boundaryMesh()[patchID]
                    );

                const edgeList& coupledPoints = cycPatch.coupledPoints();
                const labelList& meshPoints = cycPatch.meshPoints();

                label sizeby2 = cycPatch.size()/2;

//                 label ngbCycPointI = -1;
                forAll(coupledPoints, pI)
                {
                    if (localFaceID < sizeby2)
                    {
                        if ( pointI == meshPoints[coupledPoints[pI][0]] )
                        {
                            pointCyclicBndFaces[pointI] =
                                pointBndFaces
                                [
                                    meshPoints[coupledPoints[pI][1]]
                                ];
//                             ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                            break;
                        }
                    }
                    else
                    {
                        if ( pointI == meshPoints[coupledPoints[pI][1]] )
                        {
                            pointCyclicBndFaces[pointI] =
                                pointBndFaces
                                [
                                    meshPoints[coupledPoints[pI][0]]
                                ];
//                             ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                            break;
                        }
                    }
                }
            }
        }
    }

    // Find cyclicGgi boundary faces
    forAll(pointCyclicGgiBndFaces, pointI)
    {
        if (pointCyclicGgiFaces[pointI].size())
        {
            if (pointBndFaces[pointI].size())
            {
                label faceID = pointCyclicGgiFaces[pointI][0];
                label patchID = mesh().boundaryMesh().whichPatch(faceID);

                label start = mesh().boundaryMesh()[patchID].start();
                label localFaceID = faceID - start;

                const cyclicGgiPolyPatch& cycGgiPatch =
                    refCast<const cyclicGgiPolyPatch>
                    (
                        mesh().boundaryMesh()[patchID]
                    );

                const labelListList& addr
                    (
                        cycGgiPatch.master() ?
                        cycGgiPatch.patchToPatch().masterAddr()
                      : cycGgiPatch.patchToPatch().slaveAddr()
                    );

                const labelList& zoneAddr =
                    cycGgiPatch.zoneAddressing();

                label zoneLocalFaceID =
                    zoneAddr[localFaceID];

                label ngbZoneLocalFaceID =
                    addr[zoneLocalFaceID][0];

                label ngbLocalFaceID =
                    findIndex
                    (
                        cycGgiPatch.shadow().zoneAddressing(),
                        ngbZoneLocalFaceID
                    );

                if (ngbLocalFaceID == -1)
                {
                    FatalErrorIn
                    (
                        "leastSquaresVolPointInterpolation"
                        "::makePointFaces() const"
                    )
                        << "Can not find cyclic ggi patch ngb face id "
                            << mesh().boundaryMesh()[patchID].name()
                            << abort(FatalError);
                }

                const face& ngbFace =
                    cycGgiPatch.shadow().localFaces()[ngbLocalFaceID];

                const vectorField& p =
                    cycGgiPatch.shadow().localPoints();
                const labelList& meshPoints =
                    cycGgiPatch.shadow().meshPoints();

                scalar S = ngbFace.mag(p);
                scalar L = ::sqrt(S);

                label ngbPointI = -1;
                scalar minDist = GREAT;
                label minPointI = -1;
                forAll(ngbFace, pI)
                {
                    vector transCurPoint = points[pointI];

                    if (!cycGgiPatch.parallel())
                    {
                        transCurPoint =
                            transform
                            (
                                cycGgiPatch.reverseT()[0],
                                transCurPoint
                            );
                    }

                    if (cycGgiPatch.separated())
                    {
                        transCurPoint += cycGgiPatch.separationOffset();
                    }

                    scalar dist =
                        mag
                        (
                            p[ngbFace[pI]]
                          - transCurPoint
                        );

                    if (dist < 1e-2*L)
                    {
                        ngbPointI = meshPoints[ngbFace[pI]];
                        break;
                    }

                    if (dist < minDist)
                    {
                        minDist = dist;
                        minPointI = meshPoints[ngbFace[pI]];
                    }
                }

                if (ngbPointI == -1)
                {
                    FatalErrorIn
                    (
                        "leastSquaresVolPointInterpolation"
                        "::makePointFaces() const"
                    )
                        << "Can not find cyclic ggi patch ngb point "
                            << S << " " << L << " " << minDist << " "
                            << mesh().boundaryMesh()[patchID].name() << " "
                            << cycGgiPatch.master() << " "
                            << cycGgiPatch.faceCentres()[localFaceID] << " "
                            << cycGgiPatch.shadow().faceCentres()[ngbLocalFaceID] << " "
                            << points[pointI] << " " << points[minPointI] << " "
                            << transform(cycGgiPatch.reverseT()[0], points[pointI]) << " "
                            << transform(cycGgiPatch.reverseT()[0], cycGgiPatch.faceCentres()[localFaceID])
                            << abort(FatalError);
                }

                pointCyclicGgiBndFaces[pointI] = pointBndFaces[ngbPointI];
            }
        }
    }
}

void leastSquaresVolPointInterpolation::makeAxisEdges() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::"
            << "makeAxisEdges() : "
            << "constructing axis edges list"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (axisEdgesPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::"
            "makeAxisEdges() const"
        )
            << "axis edges list already exist"
            << abort(FatalError);
    }

    labelHashSet axisEdgeSet;

    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == wedgePolyPatch::typeName
        )
        {
            const wedgePolyPatch& wedge =
                refCast<const wedgePolyPatch>(mesh().boundaryMesh()[patchI]);

            const labelList& meshEdges = wedge.meshEdges();

            const labelListList& edgeFaces = mesh().edgeFaces();

            forAll(meshEdges, edgeI)
            {
                if (!wedge.isInternalEdge(edgeI))
                {
                    label curMshEdge = meshEdges[edgeI];

                    const labelList& curEdgeFaces = edgeFaces[curMshEdge];

                    if (curEdgeFaces.size() == 2)
                    {
                        label patch0 =
                            mesh().boundaryMesh().whichPatch
                            (
                                curEdgeFaces[0]
                            );

                        label patch1 =
                            mesh().boundaryMesh().whichPatch
                            (
                                curEdgeFaces[1]
                            );

                        if
                        (
                            mesh().boundaryMesh()[patch0].type()
                         == mesh().boundaryMesh()[patch1].type()
                        )
                        {
                            axisEdgeSet.insert(curMshEdge);
                        }
                    }
                }
            }

            break;
        }
    }

    axisEdgesPtr_ = new labelList(axisEdgeSet.toc());
}

void leastSquaresVolPointInterpolation::makePointAxisEdges() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::"
            << "makePointAxisEdges() : "
            << "constructing point axis edges addressing"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointAxisEdgesPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::"
            "makePointAxisEdges() const"
        )
            << "point axis edges addressing already exists"
            << abort(FatalError);
    }

    pointAxisEdgesPtr_ = new Map<labelList>();

    Map<labelList>& pointAxisEdges =
        *pointAxisEdgesPtr_;

    const edgeList& edges = mesh().edges();

    List<labelHashSet> pae(mesh().points().size());

    forAll(axisEdges(), edgeI)
    {
        label curEdge = axisEdges()[edgeI];

        pae[edges[curEdge].start()].insert(curEdge);
        pae[edges[curEdge].end()].insert(curEdge);
    }

    forAll(pae, pointI)
    {
        labelList curEdges = pae[pointI].toc();

        if (curEdges.size())
        {
            pointAxisEdges.insert(pointI, curEdges);
        }
    }

//     Info << "point-axis-edges: " << pointAxisEdges << endl;
}

// void leastSquaresVolPointInterpolation::makePointNgbProcBndFaceCentres() const
// {
//     if (debug)
//     {
//         Info<< "leastSquaresVolPointInterpolation::"
//             << "makePointNgbProcFaceCentres() : "
//             << "constructing point ngb processor face centres"
//             << endl;
//     }

//     // It is an error to attempt to recalculate
//     // if the pointer is already set
//     if (pointNgbProcBndFaceCentresPtr_)
//     {
//         FatalErrorIn
//         (
//             "leastSquaresVolPointInterpolation::"
//             "makePointNgbProcFaceCentres() const"
//         )
//             << "point ngb processor face centres already exist"
//             << abort(FatalError);
//     }

//     pointNgbProcBndFaceCentresPtr_ = new Map<vectorField>();

//     Map<vectorField>& pointNgbProcBndFaceCentres =
//         *pointNgbProcBndFaceCentresPtr_;

//     pointNgbProcBndFaceFieldData(mesh().C(), pointNgbProcBndFaceCentres);
// }

void leastSquaresVolPointInterpolation::
makeGlobalPointNgbProcBndFaceCentres() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::"
            << "makeGlobalPointNgbProcBndFaceCentres() : "
            << "constructing global point ngb processor bnd face centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalPointNgbProcBndFaceCentresPtr_)
    {
        FatalErrorIn
        (
            word("leastSquaresVolPointInterpolation::")
          + word("makeGlobalPointNgbProcBndFaceCentres() const")
        )
            << "global point ngb processor bnd face centres already exist"
                << abort(FatalError);
    }

    globalPointNgbProcBndFaceCentresPtr_ = new Map<vectorField>();

    Map<vectorField>& globalPointNgbProcBndFaceCentres =
        *globalPointNgbProcBndFaceCentresPtr_;

    globalPointNgbProcBndFaceFieldData
    (
        mesh().C(),
        globalPointNgbProcBndFaceCentres
    );
}

void leastSquaresVolPointInterpolation::
makeGlobalPointNgbProcCellCentres() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::"
            << "makeGlobalPointNgbProcCellCentres() : "
            << "constructing global point ngb processor cell centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalPointNgbProcCellCentresPtr_)
    {
        FatalErrorIn
        (
            word("leastSquaresVolPointInterpolation::")
          + word("makeGlobalPointNgbProcCellCentres() const")
        )
            << "global point ngb processor cell centres already exist"
                << abort(FatalError);
    }

    globalPointNgbProcCellCentresPtr_ = new Map<vectorField>();

    Map<vectorField>& globalPointNgbProcCellCentres =
        *globalPointNgbProcCellCentresPtr_;

    globalPointNgbProcCellFieldData(mesh().C(), globalPointNgbProcCellCentres);
}


void leastSquaresVolPointInterpolation::makeProcBndFaces() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeProcBndFaces() : "
            << "constructing list of boundary faces needed by neighbour "
            << "processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (procBndFacesPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makeProcBndFaces() const"
        )
            << "List of boundary faces needed by neighbour processors "
                << "already exist"
                << abort(FatalError);
    }

    procBndFacesPtr_ = new labelListList(Pstream::nProcs());
    labelListList& procBndFaces = *procBndFacesPtr_;
    forAll(procBndFaces, procI)
    {
        procBndFaces[procI] = labelList(0);
    }

    pointProcBndFacesPtr_ = new List<List<labelPair> >
    (
        mesh().points().size(),
        List<labelPair>(0)
    );
    List<List<labelPair> >& pointProcBndFaces = *pointProcBndFacesPtr_;

//     pointProcBndFacesPtr_ = new Map<List<labelPair> >();
//     Map<List<labelPair> >& pointProcBndFaces = *pointProcBndFacesPtr_;

    const labelListList& ptBndFaces = pointBndFaces();

//     const labelListList& pointCells = mesh().pointCells();
//     const vectorField& C = mesh().cellCentres();

    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            const labelList& bndPoints = procPatch.boundaryPoints();
            const labelList& meshPoints = procPatch.meshPoints();

            labelHashSet faceSet;
            labelHashSet pointSet;
            SLList<labelPair> pointFaceSet;

            const labelList& glPoints =
                mesh().globalData().sharedPointLabels();

            forAll(bndPoints, pointI)
            {
                label curPoint = bndPoints[pointI];
                label curMeshPoint = meshPoints[curPoint];

                bool glPoint(findIndex(glPoints, curMeshPoint) != -1);

                if (!glPoint)
                {
                    const labelList& curPointBndFaces =
                        ptBndFaces[curMeshPoint];

                    forAll(curPointBndFaces, faceI)
                    {
                        if (!faceSet.found(curPointBndFaces[faceI]))
                        {
                            faceSet.insert(curPointBndFaces[faceI]);
                        }
                        if (!pointSet.found(curPoint))
                        {
                            pointSet.insert(curPoint);
                        }
                        pointFaceSet.insert
                        (
                            labelPair(curPoint, curPointBndFaces[faceI])
                        );
                    }

//                     if (curPointBndFaces.size())
//                     {
//                         Pout << procPatch.neighbProcNo()
//                             << ", "
//                             << curMeshPoint
//                             << ", "
//                             << pointCells[curMeshPoint]
//                             << ", "
//                             << curPointBndFaces
//                             << ", "
//                             << vectorField(C, pointCells[curMeshPoint]);

//                         forAll(curPointBndFaces, fI)
//                         {
//                             Pout << ", " <<
//                                 mesh().boundaryMesh().whichPatch
//                                 (
//                                     curPointBndFaces[fI]
//                                 );
//                         }
//                         Pout << endl;
//                     }
                }
            }

            procBndFaces[procPatch.neighbProcNo()] = faceSet.toc();

            // Point face addressing
            labelList patchPoints = pointSet.toc();
            List<labelPair> patchPointsFaces(pointFaceSet);

            labelListList patchPointFaces(patchPoints.size());

            forAll(patchPoints, pointI)
            {
                labelHashSet faceSet;

                forAll(patchPointsFaces, pI)
                {
                    if
                    (
                        patchPointsFaces[pI].first()
                     == patchPoints[pointI]
                    )
                    {
                        label pointFace =
                            findIndex
                            (
                                procBndFaces[procPatch.neighbProcNo()],
                                patchPointsFaces[pI].second()
                            );
                        faceSet.insert(pointFace);
                    }
                }

                patchPointFaces[pointI] = faceSet.toc();
            }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                    // size of field
                );

                toNeighbProc << patchPoints << patchPointFaces;
            }

            labelList ngbPatchPoints;
            labelListList ngbPatchPointFaces;

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                    // size of field
                );

                fromNeighbProc >> ngbPatchPoints >> ngbPatchPointFaces;
            }

            forAll(ngbPatchPoints, pointI)
            {
                label curNgbPoint = ngbPatchPoints[pointI];

                label curLocalPoint =
                    findIndex(procPatch.neighbPoints(), curNgbPoint);
//                     procPatch.neighbPoints()[curNgbPoint];

                List<labelPair> addressing
                (
                    ngbPatchPointFaces[pointI].size(),
                    labelPair(-1, -1)
                );
                forAll(addressing, faceI)
                {
                    addressing[faceI] =
                        labelPair
                        (
                            procPatch.neighbProcNo(),
                            ngbPatchPointFaces[pointI][faceI]
                        );
                }

                pointProcBndFaces[meshPoints[curLocalPoint]] = addressing;
//                 pointProcBndFaces.insert
//                 (
//                     meshPoints[curLocalPoint],
//                     addressing
//                 );
            }
        }
    }
}


void leastSquaresVolPointInterpolation::makeProcBndFaceCentres() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeProcBndFaceCentres() : "
            << "constructing centres of boundary faces from ngb processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (procBndFaceCentresPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makeProcBndFaceCentres() const"
        )
            << "Centres of faces from ngb processors already exist"
                << abort(FatalError);
    }

    procBndFaceCentresPtr_ =
        new FieldField<Field, vector>
        (
            Pstream::nProcs()
        );
    FieldField<Field, vector>& procBndFaceCentres = *procBndFaceCentresPtr_;

    const vectorField& Cf = mesh().faceCentres();

    if (Pstream::parRun())
    {
        // Send centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres(Cf, procBndFaces()[procI]);

                {
                    OPstream toNeighbProc
                    (
                        Pstream::blocking,
                        procI
                        // size of field
                    );

                    toNeighbProc << procCentres;
                }
            }
        }

        // Receive centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres;

                {
                    IPstream fromNeighbProc
                    (
                        Pstream::blocking,
                        procI
                        // size of field
                    );

                    fromNeighbProc >> procCentres;
                }

                procBndFaceCentres.set(procI, new vectorField(procCentres));
            }
            else
            {
                procBndFaceCentres.set(procI, new vectorField(0));
            }
        }
    }
    else
    {
        forAll (procBndFaceCentres, procI)
        {
            procBndFaceCentres.set(procI, new vectorField(0));
        }
    }
}


void leastSquaresVolPointInterpolation::makeProcCells() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeProcCells() : "
            << "constructing list of cells needed by neighbour processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (procCellsPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makeProcCells() const"
        )
            << "List of cells needed by neighbour processors already exist"
                << abort(FatalError);
    }

    procCellsPtr_ = new labelListList(Pstream::nProcs());
    labelListList& procCells = *procCellsPtr_;
    forAll(procCells, procI)
    {
        procCells[procI] = labelList(0);
    }

    pointProcCellsPtr_ = new Map<List<labelPair> >();
    Map<List<labelPair> >& pointProcCells = *pointProcCellsPtr_;

    const labelListList& pointCells = mesh().pointCells();

    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            const labelList& meshPoints = procPatch.meshPoints();

            const unallocLabelList& patchCells = procPatch.faceCells();

            labelHashSet cellSet;
            labelHashSet pointSet;
            SLList<labelPair> pointCellSet;

            const labelList& glPoints =
                mesh().globalData().sharedPointLabels();

            forAll(meshPoints, pointI)
            {
                label curPoint = meshPoints[pointI];

                bool glPoint(findIndex(glPoints, curPoint) != -1);

                labelHashSet localCellSet;

                if (!glPoint)
                {
                    const labelList& curCells = pointCells[curPoint];

                    forAll(curCells, cellI)
                    {
                        if (findIndex(patchCells, curCells[cellI]) == -1)
                        {
                            if (!cellSet.found(curCells[cellI]))
                            {
                                cellSet.insert(curCells[cellI]);
                            }

                            if (!localCellSet.found(curCells[cellI]))
                            {
                                localCellSet.insert(curCells[cellI]);
                            }

                            if (!pointSet.found(pointI))
                            {
                                pointSet.insert(pointI);
                            }

                            pointCellSet.insert
                            (
                                labelPair(pointI, curCells[cellI])
                            );
                        }
                    }
                }
            }

            procCells[procPatch.neighbProcNo()] = cellSet.toc();

            // Point cell addressing
            labelList patchPoints = pointSet.toc();
            List<labelPair> patchPointsCells(pointCellSet);

            labelListList patchPointCells(patchPoints.size());

            label nCells = 0;

            forAll(patchPointCells, pointI)
            {
                labelHashSet cellSet;

                forAll(patchPointsCells, pI)
                {
                    if
                    (
                        patchPointsCells[pI].first()
                     == patchPoints[pointI]
                    )
                    {
                        label pointCell =
                            findIndex
                            (
                                procCells[procPatch.neighbProcNo()],
                                patchPointsCells[pI].second()
                            );
                        cellSet.insert(pointCell);
                    }
                }

                patchPointCells[pointI] = cellSet.toc();

                nCells += patchPointCells[pointI].size();
            }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                    // size of field
//                   + 2*patchPoints.size()*sizeof(label)
//                   + nCells*sizeof(label)
                );

                toNeighbProc << patchPoints << patchPointCells;
            }
        }
    }


    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            const labelList& meshPoints = procPatch.meshPoints();

            labelList ngbPatchPoints;
            labelListList ngbPatchPointCells;

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                    // size of field
                );

                fromNeighbProc >> ngbPatchPoints >> ngbPatchPointCells;
            }

            forAll(ngbPatchPoints, pointI)
            {
                label curNgbPoint = ngbPatchPoints[pointI];

                label curLocalPoint =
                    findIndex(procPatch.neighbPoints(), curNgbPoint);
//                     procPatch.neighbPoints()[curNgbPoint];

                List<labelPair> addressing
                (
                    ngbPatchPointCells[pointI].size(),
                    labelPair(-1, -1)
                );

                forAll(addressing, cellI)
                {
                    addressing[cellI] =
                        labelPair
                        (
                            procPatch.neighbProcNo(),
                            ngbPatchPointCells[pointI][cellI]
                        );
                }

                pointProcCells.insert
                (
                    meshPoints[curLocalPoint],
                    addressing
                );
            }
        }
    }
}


void leastSquaresVolPointInterpolation::makeProcCellCentres() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeProcCellCentres() : "
            << "constructing centres of cells from ngb processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (procCellCentresPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makeProcCellCentres() const"
        )
            << "Centres of cells from ngb processors already exist"
                << abort(FatalError);
    }

    procCellCentresPtr_ =
        new FieldField<Field, vector>
        (
            Pstream::nProcs()
        );
    FieldField<Field, vector>& procCellCentres = *procCellCentresPtr_;

//     const vectorField& CI = mesh().C().internalField();
    const vectorField& CI = mesh().cellCentres();

    if (Pstream::parRun())
    {
        // Send centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres(CI, procCells()[procI]);

                {
                    OPstream toNeighbProc
                    (
                        Pstream::blocking,
                        procI
                        // size of field
                    );

                    toNeighbProc << procCentres;
                }
            }
        }

        // Receive centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres;

                {
                    IPstream fromNeighbProc
                    (
                        Pstream::blocking,
                        procI
                        // size of field
                    );

                    fromNeighbProc >> procCentres;
                }

                procCellCentres.set(procI, new vectorField(procCentres));
            }
            else
            {
                procCellCentres.set(procI, new vectorField(0));
            }
        }
    }
    else
    {
        forAll (procCellCentres, procI)
        {
            procCellCentres.set(procI, new vectorField(0));
        }
    }
}


void leastSquaresVolPointInterpolation::makeWeights() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeWeights() : "
            << "constructing weights"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (weightsPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makeWeights() const"
        )
            << "weights already exist"
            << abort(FatalError);
    }

    weightsPtr_ =
        new FieldField<Field, scalar>(mesh().points().size());
    FieldField<Field, scalar>& weights = *weightsPtr_;

    const vectorField& p = mesh().points();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Cf = mesh().faceCentres();

    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptCyclicGgiFaces = pointCyclicGgiFaces();
    const labelListList& ptCyclicGgiBndFaces = pointCyclicGgiBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

//     const Map<vectorField>& ptNgbProcBndFaceCentres =
//         pointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcCellCentres =
        globalPointNgbProcCellCentres();

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, vector>& procCentres = procCellCentres();

    const FieldField<Field, vector>& procBndFaceCent = procBndFaceCentres();

    forAll(weights, pointI)
    {
        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpCyclicGgiFaces = ptCyclicGgiFaces[pointI];
        const labelList& interpCyclicGgiBndFaces = ptCyclicGgiBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

//         vectorField interpNgbProcBndFaceCentres(0);

//         // Boundary faces from neighbour processors
//         if (ptNgbProcBndFaceCentres.found(pointI))
//         {
//              interpNgbProcBndFaceCentres =
//                 ptNgbProcBndFaceCentres[pointI];
//         }

        vectorField glInterpNgbProcBndFaceCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcBndFaceCentres.found(pointI))
        {
             glInterpNgbProcBndFaceCentres =
                gPtNgbProcBndFaceCentres[pointI];
        }

        vectorField glInterpNgbProcCellCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcCellCentres.found(pointI))
        {
             glInterpNgbProcCellCentres =
                gPtNgbProcCellCentres[pointI];
        }

        vectorField interpNgbProcCellCentres(0);

        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellCentres.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellCentres[cI] =
                    procCentres
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        vectorField interpNgbProcBndFaceCentres(0);

        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceCentres.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceCentres[fI] =
                    procBndFaceCent
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        vectorField allPoints
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpCyclicGgiFaces.size()
          + interpCyclicGgiBndFaces.size()
          + interpProcFaces.size()
//           + interpNgbProcBndFaceCentres.size()
          + glInterpNgbProcBndFaceCentres.size()
          + glInterpNgbProcCellCentres.size()
          + interpNgbProcCellCentres.size()
          + interpNgbProcBndFaceCentres.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Boundary faces
        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];

            allPoints[pointID++] = Cf[faceID];
        }

        // Cyclic boundary faces
        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if (cycPatch.parallel())
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaOwn - deltaNgb;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaNgb - deltaOwn;
                    delta *= -1;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = // dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaOwn
                      - transform(cycPatch.forwardT()[0], deltaNgb);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaNgb
                      - transform(cycPatch.forwardT()[0], deltaOwn);

                    delta = -transform(cycPatch.reverseT()[0], delta);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[cycPatchID]);

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];

            const edgeList& coupledPoints = cycPolyPatch.coupledPoints();
            const labelList& meshPoints = cycPolyPatch.meshPoints();

            label sizeby2 = cycPolyPatch.size()/2;

            label ngbCycPointI = -1;
            forAll(coupledPoints, pI)
            {
                if (cycLocalFaceID < sizeby2)
                {
                    if ( pointI == meshPoints[coupledPoints[pI][0]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                        break;
                    }
                }
                else
                {
                    if ( pointI == meshPoints[coupledPoints[pI][1]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                        break;
                    }
                }
            }

            if (cycPatch.parallel())
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];

                allPoints[pointID++] = p[pointI] + deltaNgb;
            }
            else
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];


                if (cycLocalFaceID < sizeby2)
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.forwardT()[0], deltaNgb);
                }
                else
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.reverseT()[0], deltaNgb);
                }
            }
        }

        // Cyclic ggi boundary faces
        for (label i=0; i<interpCyclicGgiFaces.size(); i++)
        {
            label faceID = interpCyclicGgiFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[patchID]
                );

            const vectorField& rfcc =
                cycGgiPatch.reconFaceCellCentres();

            allPoints[pointID++] = rfcc[localFaceID];
        }

        if (interpCyclicGgiBndFaces.size())
        {
            label cycGgiFaceID = interpCyclicGgiFaces[0];
            label cycGgiPatchID =
                mesh().boundaryMesh().whichPatch(cycGgiFaceID);
            label cycGgiStart = mesh().boundaryMesh()[cycGgiPatchID].start();
            label cycGgiLocalFaceID = cycGgiFaceID - cycGgiStart;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[cycGgiPatchID]
                );

            const labelListList& addr
                (
                    cycGgiPatch.master() ?
                    cycGgiPatch.patchToPatch().masterAddr()
                  : cycGgiPatch.patchToPatch().slaveAddr()
                );

            const labelList& zoneAddr =
                cycGgiPatch.zoneAddressing();

            label zoneLocalFaceID = zoneAddr[cycGgiLocalFaceID];

            label ngbZoneLocalFaceID = addr[zoneLocalFaceID][0];

            label ngbLocalFaceID =
                findIndex
                (
                    cycGgiPatch.shadow().zoneAddressing(),
                    ngbZoneLocalFaceID
                );
            if (ngbLocalFaceID == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makeWeights() const"
                )
                    << "Can not find cyclic ggi patch ngb face id"
                        << abort(FatalError);
            }

            const face& ngbFace =
                cycGgiPatch.shadow().localFaces()[ngbLocalFaceID];

            const vectorField& lp =
                cycGgiPatch.shadow().localPoints();
            const labelList& meshPoints =
                cycGgiPatch.shadow().meshPoints();

            scalar S = ngbFace.mag(lp);
            scalar L = ::sqrt(S);

            label ngbPointI = -1;
            forAll(ngbFace, pI)
            {
                vector transCurPoint = p[pointI];

                if (!cycGgiPatch.parallel())
                {
                    transCurPoint =
                        transform
                        (
                            cycGgiPatch.reverseT()[0],
                            transCurPoint
                        );
                }

                if (cycGgiPatch.separated())
                {
                    transCurPoint += cycGgiPatch.separationOffset();
                }

                scalar dist =
                    mag
                    (
                        lp[ngbFace[pI]]
                      - transCurPoint
                    );

                if (dist < 1e-2*L)
                {
                    ngbPointI = meshPoints[ngbFace[pI]];
                    break;
                }
            }

            if (ngbPointI == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makeWeights() const"
                )
                    << "Can not find cyclic ggi patch ngb point"
                        << abort(FatalError);
            }


            for (label i=0; i<interpCyclicGgiBndFaces.size(); i++)
            {
                label faceID = interpCyclicGgiBndFaces[i];

                if (cycGgiPatch.parallel())
                {
                    vector deltaNgb = //dni
                        Cf[faceID] - p[ngbPointI];

                    allPoints[pointID++] = p[pointI] + deltaNgb;
                }
                else
                {
                    vector deltaNgb = //dni
                        Cf[faceID] - p[ngbPointI];

                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycGgiPatch.forwardT()[0], deltaNgb);
                }
            }
        }

        // Processor boundary faces
        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            allPoints[pointID++] =
                mesh().C().boundaryField()[patchID][localFaceID];
        }

//         // Boundary faces from neighbour processors
//         for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
//         {
//             allPoints[pointID++] =
//                 interpNgbProcBndFaceCentres[i];
//         }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcBndFaceCentres[i];
        }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcCellCentres[i];
        }

        // Cells from neighbour processors
        for (label i=0; i<interpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcCellCentres[i];
        }

        // Boundary faces from neighbour processors
        for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcBndFaceCentres[i];
        }

        vectorField allMirrorPoints(0);
        if (mag(mirrorPlaneTransformation()[pointI].first())>SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const vector& n = mirrorPlaneTransformation()[pointI].first();

            allMirrorPoints.setSize(allPoints.size());

            forAll(allPoints, pI)
            {
                allMirrorPoints[pI] =
                    p[pointI] + transform(I-2*n*n, (allPoints[pI]-p[pointI]));
            }
        }

        // Weights
        scalarField W(allPoints.size() + allMirrorPoints.size(), 1.0);

//         label pI = 0;
//         for (label i=0; i<allPoints.size(); i++)
//         {
//             scalar curR =  mag(allPoints[i] - p[pointI]);
//             W[pI++] = 1.0/(sqr(curR) + VSMALL);
//         }
//         for (label i=0; i<allMirrorPoints.size(); i++)
//         {
//             scalar curR =  mag(allMirrorPoints[i] - p[pointI]);
//             W[pI++] = 1.0/(sqr(curR) + VSMALL);
//         }

        weights.set(pointI, new scalarField(W));
    }
}


void leastSquaresVolPointInterpolation::makeOrigins() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeOrigin() : "
            << "constructing local origins"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (originsPtr_ || refLPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::makeOrigin() const"
        )
            << "local origins (or refLs) already exist"
            << abort(FatalError);
    }

    originsPtr_ = new vectorField(mesh().points().size(), vector::zero);
    vectorField& origins = *originsPtr_;

    refLPtr_ = new scalarField(mesh().points().size(), 0);
    scalarField& refL = *refLPtr_;


    const FieldField<Field, scalar>& w = weights();

    const vectorField& p = mesh().points();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Cf = mesh().faceCentres();

    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptCyclicGgiFaces = pointCyclicGgiFaces();
    const labelListList& ptCyclicGgiBndFaces = pointCyclicGgiBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

//     const Map<vectorField> ptNgbProcBndFaceCentres =
//         pointNgbProcBndFaceCentres();

    const Map<vectorField> gPtNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentres();

    const Map<vectorField> gPtNgbProcCellCentres =
        globalPointNgbProcCellCentres();

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, vector>& procCentres = procCellCentres();
    const FieldField<Field, vector>& procBndFaceCent = procBndFaceCentres();

    forAll(origins, pointI)
    {
        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpCyclicGgiFaces = ptCyclicGgiFaces[pointI];
        const labelList& interpCyclicGgiBndFaces = ptCyclicGgiBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

//         vectorField interpNgbProcBndFaceCentres(0);

//         // Boundary faces from neighbour processors
//         if (ptNgbProcBndFaceCentres.found(pointI))
//         {
//              interpNgbProcBndFaceCentres =
//                 ptNgbProcBndFaceCentres[pointI];
//         }

        vectorField glInterpNgbProcBndFaceCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcBndFaceCentres.found(pointI))
        {
             glInterpNgbProcBndFaceCentres =
                gPtNgbProcBndFaceCentres[pointI];
        }

        vectorField glInterpNgbProcCellCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcCellCentres.found(pointI))
        {
             glInterpNgbProcCellCentres =
                gPtNgbProcCellCentres[pointI];
        }

        vectorField interpNgbProcCellCentres(0);

        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellCentres.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellCentres[cI] =
                    procCentres
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        vectorField interpNgbProcBndFaceCentres(0);

        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceCentres.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceCentres[fI] =
                    procBndFaceCent
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        vectorField allPoints
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpCyclicGgiFaces.size()
          + interpCyclicGgiBndFaces.size()
          + interpProcFaces.size()
//           + interpNgbProcBndFaceCentres.size()
          + glInterpNgbProcBndFaceCentres.size()
          + glInterpNgbProcCellCentres.size()
          + interpNgbProcCellCentres.size()
          + interpNgbProcBndFaceCentres.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Boundary faces
        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];

            allPoints[pointID++] = Cf[faceID];
        }

        // Cyclic boundary faces
        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if (cycPatch.parallel())
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaOwn - deltaNgb;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaNgb - deltaOwn;
                    delta *= -1;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = // dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaOwn
                      - transform(cycPatch.forwardT()[0], deltaNgb);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "own: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaNgb
                      - transform(cycPatch.forwardT()[0], deltaOwn);

                    delta = -transform(cycPatch.reverseT()[0], delta);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "ngb: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
            }

//             label faceID = interpCyclicFaces[i];
//             label patchID = mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             const unallocLabelList& faceCells =
//                 mesh().boundary()[patchID].faceCells();

//             label sizeby2 = faceCells.size()/2;

//             if (localFaceID < sizeby2)
//             {
//                 vector delta =
//                     C[faceCells[localFaceID + sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID + sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
//             else
//             {
//                 vector delta =
//                     C[faceCells[localFaceID - sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID - sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[cycPatchID]);

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];

            const edgeList& coupledPoints = cycPolyPatch.coupledPoints();
            const labelList& meshPoints = cycPolyPatch.meshPoints();

            label sizeby2 = cycPolyPatch.size()/2;

            label ngbCycPointI = -1;
            forAll(coupledPoints, pI)
            {
                if (cycLocalFaceID < sizeby2)
                {
                    if ( pointI == meshPoints[coupledPoints[pI][0]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                        break;
                    }
                }
                else
                {
                    if ( pointI == meshPoints[coupledPoints[pI][1]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                        break;
                    }
                }
            }

            if (cycPatch.parallel())
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];

                allPoints[pointID++] = p[pointI] + deltaNgb;
            }
            else
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];


                if (cycLocalFaceID < sizeby2)
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.forwardT()[0], deltaNgb);
                }
                else
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.reverseT()[0], deltaNgb);
                }
            }
        }

        // Cyclic ggi boundary faces
        for (label i=0; i<interpCyclicGgiFaces.size(); i++)
        {
            label faceID = interpCyclicGgiFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[patchID]
                );

            const vectorField& rfcc =
                cycGgiPatch.reconFaceCellCentres();

            allPoints[pointID++] = rfcc[localFaceID];
        }

        if (interpCyclicGgiBndFaces.size())
        {
            label cycGgiFaceID = interpCyclicGgiFaces[0];
            label cycGgiPatchID =
                mesh().boundaryMesh().whichPatch(cycGgiFaceID);
            label cycGgiStart = mesh().boundaryMesh()[cycGgiPatchID].start();
            label cycGgiLocalFaceID = cycGgiFaceID - cycGgiStart;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[cycGgiPatchID]
                );

            const labelListList& addr
                (
                    cycGgiPatch.master() ?
                    cycGgiPatch.patchToPatch().masterAddr()
                  : cycGgiPatch.patchToPatch().slaveAddr()
                );

            const labelList& zoneAddr =
                cycGgiPatch.zoneAddressing();

            label zoneLocalFaceID = zoneAddr[cycGgiLocalFaceID];

            label ngbZoneLocalFaceID = addr[zoneLocalFaceID][0];

            label ngbLocalFaceID =
                findIndex
                (
                    cycGgiPatch.shadow().zoneAddressing(),
                    ngbZoneLocalFaceID
                );
            if (ngbLocalFaceID == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makeOrigins() const"
                )
                    << "Can not find cyclic ggi patch ngb face id"
                        << abort(FatalError);
            }

            const face& ngbFace =
                cycGgiPatch.shadow().localFaces()[ngbLocalFaceID];

            const vectorField& lp =
                cycGgiPatch.shadow().localPoints();
            const labelList& meshPoints =
                cycGgiPatch.shadow().meshPoints();

            scalar S = ngbFace.mag(lp);
            scalar L = ::sqrt(S);

            label ngbPointI = -1;
            forAll(ngbFace, pI)
            {
                vector transCurPoint = p[pointI];

                if (!cycGgiPatch.parallel())
                {
                    transCurPoint =
                        transform
                        (
                            cycGgiPatch.reverseT()[0],
                            transCurPoint
                        );
                }

                if (cycGgiPatch.separated())
                {
                    transCurPoint += cycGgiPatch.separationOffset();
                }

                scalar dist =
                    mag
                    (
                        lp[ngbFace[pI]]
                      - transCurPoint
                    );

//                 scalar dist =
//                     mag
//                     (
//                         lp[ngbFace[pI]]
//                       - transform
//                         (
//                             cycGgiPatch.reverseT()[0],
//                             p[pointI]
//                         )
//                     );
                if (dist < 1e-2*L)
                {
                    ngbPointI = meshPoints[ngbFace[pI]];
                    break;
                }
            }
            if (ngbPointI == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makeOrigins() const"
                )
                    << "Can not find cyclic ggi patch ngb point"
                        << abort(FatalError);
            }


            for (label i=0; i<interpCyclicGgiBndFaces.size(); i++)
            {
                label faceID = interpCyclicGgiBndFaces[i];

                if (cycGgiPatch.parallel())
                {
                    vector deltaNgb = //dni
                        Cf[faceID] - p[ngbPointI];

                    allPoints[pointID++] = p[pointI] + deltaNgb;
                }
                else
                {
                    vector deltaNgb = //dni
                        Cf[faceID] - p[ngbPointI];

                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycGgiPatch.forwardT()[0], deltaNgb);
                }
            }
        }

        // Processor boundary faces
        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            allPoints[pointID++] =
                mesh().C().boundaryField()[patchID][localFaceID];
        }

//         // Boundary faces from neighbour processors
//         for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
//         {
//             allPoints[pointID++] =
//                 interpNgbProcBndFaceCentres[i];
//         }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcBndFaceCentres[i];
        }

        // Global point cells from neighbour processors
        for (label i=0; i<glInterpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcCellCentres[i];
        }

        // Cells from neighbour processors
        for (label i=0; i<interpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcCellCentres[i];
        }

        // Boundary faces from neighbour processors
        for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcBndFaceCentres[i];
        }

        vectorField allMirrorPoints(0);
        if (mag(mirrorPlaneTransformation()[pointI].first())>SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const vector& n = mirrorPlaneTransformation()[pointI].first();

            allMirrorPoints.setSize(allPoints.size());

            forAll(allPoints, pI)
            {
                allMirrorPoints[pI] =
                    p[pointI] + transform(I-2*n*n, (allPoints[pI]-p[pointI]));
            }
        }


        const scalarField& W = w[pointI];

        label pI = 0;
        for (label i=0; i<allPoints.size(); i++)
        {
            origins[pointI] += sqr(W[pI++])*allPoints[i];
        }
        for (label i=0; i<allMirrorPoints.size(); i++)
        {
            origins[pointI] += sqr(W[pI++])*allMirrorPoints[i];
        }

        origins[pointI] /= sum(sqr(W));

        boundBox bb(allPoints, false);
        refL[pointI] = mag(bb.max() - bb.min())/2;
    }
}


void leastSquaresVolPointInterpolation::makeInvLsMatrices() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::makeInvLsMatrices() : "
            << "making least squares linear interpolation matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invLsMatrices_.size() != 0)
    {
        FatalErrorIn("leastSquaresVolPointInterpolation::makeInvLsMatrices()")
            << "least square linear inerpolation matrices already exist"
            << abort(FatalError);
    }

    invLsMatrices_.setSize(mesh().points().size());

    const vectorField& p = mesh().points();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Cf = mesh().faceCentres();

    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptCyclicGgiFaces = pointCyclicGgiFaces();
    const labelListList& ptCyclicGgiBndFaces = pointCyclicGgiBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

//     const Map<vectorField>& ptNgbProcBndFaceCentres =
//         pointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcCellCentres =
        globalPointNgbProcCellCentres();

    const Map<List<labelPair> >& ptProcCells = pointProcCells();
    const FieldField<Field, vector>& procCentres = procCellCentres();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
    const FieldField<Field, vector>& procBndFaceCent = procBndFaceCentres();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();

    const scalarField& L = refL();

    label nCoeffs = 3;

    scalarField D(invLsMatrices_.size(), 0);

    forAll(invLsMatrices_, pointI)
    {
        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpCyclicGgiFaces = ptCyclicGgiFaces[pointI];
        const labelList& interpCyclicGgiBndFaces = ptCyclicGgiBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

//         vectorField interpNgbProcBndFaceCentres(0);// = vectorField::null();

//         // Boundary faces from neighbour processors
//         if (ptNgbProcBndFaceCentres.found(pointI))
//         {
//              interpNgbProcBndFaceCentres =
//                 ptNgbProcBndFaceCentres[pointI];
//         }

        vectorField glInterpNgbProcBndFaceCentres(0);

        // Boundar faces from neighbour processors
        if (gPtNgbProcBndFaceCentres.found(pointI))
        {
             glInterpNgbProcBndFaceCentres =
                gPtNgbProcBndFaceCentres[pointI];
        }

        vectorField glInterpNgbProcCellCentres(0);

        // Cells from neighbour processors
        if (gPtNgbProcCellCentres.found(pointI))
        {
             glInterpNgbProcCellCentres =
                gPtNgbProcCellCentres[pointI];
        }

        vectorField interpNgbProcCellCentres(0);

        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellCentres.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellCentres[cI] =
                    procCentres
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        vectorField interpNgbProcBndFaceCentres(0);

        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceCentres.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceCentres[fI] =
                    procBndFaceCent
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        vectorField allPoints
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpCyclicGgiFaces.size()
          + interpCyclicGgiBndFaces.size()
          + interpProcFaces.size()
//           + interpNgbProcBndFaceCentres.size()
          + glInterpNgbProcBndFaceCentres.size()
          + glInterpNgbProcCellCentres.size()
          + interpNgbProcCellCentres.size()
          + interpNgbProcBndFaceCentres.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Boundary faces
        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];

            allPoints[pointID++] = Cf[faceID];
        }

        // Cyclic boundary faces
        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if (cycPatch.parallel())
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaOwn - deltaNgb;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaNgb - deltaOwn;
                    delta *= -1;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = // dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaOwn
                      - transform(cycPatch.forwardT()[0], deltaNgb);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "own: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaNgb
                      - transform(cycPatch.forwardT()[0], deltaOwn);

                    delta = -transform(cycPatch.reverseT()[0], delta);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "ngb: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
            }

//             label faceID = interpCyclicFaces[i];
//             label patchID = mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             const unallocLabelList& faceCells =
//                 mesh().boundary()[patchID].faceCells();

//             label sizeby2 = faceCells.size()/2;

//             if (localFaceID < sizeby2)
//             {
//                 vector delta =
//                     C[faceCells[localFaceID + sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID + sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
//             else
//             {
//                 vector delta =
//                     C[faceCells[localFaceID - sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID - sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[cycPatchID]);

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];

            const edgeList& coupledPoints = cycPolyPatch.coupledPoints();
            const labelList& meshPoints = cycPolyPatch.meshPoints();

            label sizeby2 = cycPolyPatch.size()/2;

            label ngbCycPointI = -1;
            forAll(coupledPoints, pI)
            {
                if (cycLocalFaceID < sizeby2)
                {
                    if ( pointI == meshPoints[coupledPoints[pI][0]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                        break;
                    }
                }
                else
                {
                    if ( pointI == meshPoints[coupledPoints[pI][1]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                        break;
                    }
                }
            }

            if (cycPatch.parallel())
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];

                allPoints[pointID++] = p[pointI] + deltaNgb;
            }
            else
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];


                if (cycLocalFaceID < sizeby2)
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.forwardT()[0], deltaNgb);
                }
                else
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.reverseT()[0], deltaNgb);
                }
            }
        }

        // Cyclic ggi boundary faces
        for (label i=0; i<interpCyclicGgiFaces.size(); i++)
        {
            label faceID = interpCyclicGgiFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[patchID]
                );

            const vectorField& rfcc =
                cycGgiPatch.reconFaceCellCentres();

            allPoints[pointID++] = rfcc[localFaceID];
        }

        if (interpCyclicGgiBndFaces.size())
        {
            label cycGgiFaceID = interpCyclicGgiFaces[0];
            label cycGgiPatchID =
                mesh().boundaryMesh().whichPatch(cycGgiFaceID);
            label cycGgiStart = mesh().boundaryMesh()[cycGgiPatchID].start();
            label cycGgiLocalFaceID = cycGgiFaceID - cycGgiStart;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[cycGgiPatchID]
                );

            const labelListList& addr
                (
                    cycGgiPatch.master() ?
                    cycGgiPatch.patchToPatch().masterAddr()
                  : cycGgiPatch.patchToPatch().slaveAddr()
                );

            const labelList& zoneAddr =
                cycGgiPatch.zoneAddressing();

            label zoneLocalFaceID = zoneAddr[cycGgiLocalFaceID];

            label ngbZoneLocalFaceID = addr[zoneLocalFaceID][0];

            label ngbLocalFaceID =
                findIndex
                (
                    cycGgiPatch.shadow().zoneAddressing(),
                    ngbZoneLocalFaceID
                );
            if (ngbLocalFaceID == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makeInvLsMatrices() const"
                )
                    << "Can not find cyclic ggi patch ngb face id"
                        << abort(FatalError);
            }

            const face& ngbFace =
                cycGgiPatch.shadow().localFaces()[ngbLocalFaceID];

            const vectorField& lp =
                cycGgiPatch.shadow().localPoints();
            const labelList& meshPoints =
                cycGgiPatch.shadow().meshPoints();

            scalar S = ngbFace.mag(lp);
            scalar L = ::sqrt(S);

            label ngbPointI = -1;
            forAll(ngbFace, pI)
            {
                vector transCurPoint = p[pointI];

                if (!cycGgiPatch.parallel())
                {
                    transCurPoint =
                        transform
                        (
                            cycGgiPatch.reverseT()[0],
                            transCurPoint
                        );
                }

                if (cycGgiPatch.separated())
                {
                    transCurPoint += cycGgiPatch.separationOffset();
                }

                scalar dist =
                    mag
                    (
                        lp[ngbFace[pI]]
                      - transCurPoint
                    );

//                 scalar dist =
//                     mag
//                     (
//                         lp[ngbFace[pI]]
//                       - transform
//                         (
//                             cycGgiPatch.reverseT()[0],
//                             p[pointI]
//                         )
//                     );

                if (dist < 1e-2*L)
                {
                    ngbPointI = meshPoints[ngbFace[pI]];
                    break;
                }
            }
            if (ngbPointI == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makeInvLsMatrices() const"
                )
                    << "Can not find cyclic ggi patch ngb point"
                        << abort(FatalError);
            }

            for (label i=0; i<interpCyclicGgiBndFaces.size(); i++)
            {
                label faceID = interpCyclicGgiBndFaces[i];

                if (cycGgiPatch.parallel())
                {
                    vector deltaNgb = //dni
                        Cf[faceID] - p[ngbPointI];

                    allPoints[pointID++] = p[pointI] + deltaNgb;
                }
                else
                {
                    vector deltaNgb = //dni
                        Cf[faceID] - p[ngbPointI];

                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycGgiPatch.forwardT()[0], deltaNgb);
                }
            }
        }

        // Processor boundary faces
        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            allPoints[pointID++] =
                mesh().C().boundaryField()[patchID][localFaceID];
        }

//         // Boundary faces from neighbour processors
//         for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
//         {
//             allPoints[pointID++] =
//                 interpNgbProcBndFaceCentres[i];
//         }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcBndFaceCentres[i];
        }

        // Global point cells from neighbour processors
        for (label i=0; i<glInterpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcCellCentres[i];
        }

        // Cells from neighbour processors
        for (label i=0; i<interpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcCellCentres[i];
        }

        // Bnd faces from neighbour processors
        for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcBndFaceCentres[i];
        }

        vectorField allMirrorPoints(0);
        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const vector& n = mirrorPlaneTransformation()[pointI].first();

            allMirrorPoints.setSize(allPoints.size());

            forAll(allPoints, pI)
            {
                allMirrorPoints[pI] =
                    p[pointI] + transform(I-2*n*n, (allPoints[pI]-p[pointI]));
            }
        }

        if (allPoints.size() + allMirrorPoints.size() < nCoeffs)
        {
            Pout << pointI << ", "
                << interpCells.size() << ", "
                << interpBndFaces.size() << ", "
                << interpCyclicFaces.size() << ", "
                << interpCyclicBndFaces.size() << ", "
                << interpProcFaces.size() << ", "
                << interpNgbProcCellCentres.size() << endl;

            FatalErrorIn
            (
                "leastSquaresVolPointInterpolation::makeInvInvMatrices()"
            )
                << "allPoints.size() < " << nCoeffs << " : "
                    << allPoints.size() + allMirrorPoints.size()
                    << abort(FatalError);
        }

//         // Weights
//         scalarField W(allPoints.size(), 1.0);
//         scalar sumW = 0;
//         for (label i=0; i<allPoints.size(); i++)
//         {
//             scalar curR =  mag(allPoints[i] - p[pointI]);
//             W[i] = 1.0/(sqr(curR) + VSMALL);
//             sumW += W[i];
//         }
//         W /= sumW;

        const scalarField& W = w[pointI];

        invLsMatrices_.set
        (
            pointI,
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size() + allMirrorPoints.size(),
                0.0
            )
        );
        scalarRectangularMatrix& curMatrix = invLsMatrices_[pointI];

        scalarRectangularMatrix M
        (
            allPoints.size() + allMirrorPoints.size(),
            nCoeffs,
            0.0
        );

        label pI = 0;
        for (label i=0; i<allPoints.size(); i++)
        {
            scalar X = (allPoints[i].x() - o[pointI].x())/L[pointI];
            scalar Y = (allPoints[i].y() - o[pointI].y())/L[pointI];
            scalar Z = (allPoints[i].z() - o[pointI].z())/L[pointI];

            M[pI][0] = X;
            M[pI][1] = Y;
            M[pI][2] = Z;
            pI++;
        }
        for (label i=0; i<allMirrorPoints.size(); i++)
        {
            scalar X = (allMirrorPoints[i].x() - o[pointI].x())/L[pointI];
            scalar Y = (allMirrorPoints[i].y() - o[pointI].y())/L[pointI];
            scalar Z = (allMirrorPoints[i].z() - o[pointI].z())/L[pointI];

            M[pI][0] = X;
            M[pI][1] = Y;
            M[pI][2] = Z;
            pI++;
        }

        // Applying weights
        for (label i=0; i<M.n(); i++)
        {
            for (label j=0; j<M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

//         SVD svd(M, SMALL);

//         for (label i=0; i<svd.VSinvUt().n(); i++)
//         {
//             for (label j=0; j<svd.VSinvUt().m(); j++)
//             {
//                 curMatrix[i][j] = svd.VSinvUt()[i][j]*W[j];
//             }
//         }

//         scalarSquareMatrix lsM(nCoeffs, 0.0);
        tensor lsM = tensor::zero;

        for (label i=0; i<3; i++)
        {
            for (label j=0; j<3; j++)
            {
                for (label k=0; k<M.n(); k++)
                {
//                     lsM[i][j] += M[k][i]*M[k][j];
                    lsM(i,j) += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate matrix norm
        scalar maxRowSum = 0.0;
        for (label i=0; i<3; i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<3; j++)
            {
//                 curRowSum += lsM[i][j];
                curRowSum += lsM(i,j);
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }

        // Calculate inverse
//         scalarSquareMatrix invLsM = lsM.LUinvert();

        D[pointI] = det(lsM);

        if (mag(D[pointI]) > SMALL)
        {
            tensor invLsM = hinv(lsM); // hinv(lsM)

            for (label i=0; i<3; i++)
            {
                for (label j=0; j<M.n(); j++)
                {
                    for (label k=0; k<3; k++)
                    {
//                     curMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
                        curMatrix[i][j] += invLsM(i,k)*M[j][k]*W[j];
                    }
                }
            }
        }
        else
        {
            Pout << "Det: " << D[pointI] << endl;

            for (label i=0; i<3; i++)
            {
                for (label j=0; j<M.n(); j++)
                {
                    curMatrix[i][j] = 0;
                }
            }
        }
    }

    // Pout << "det, min: " << min(D) << ", avg: " << average(D)
    //     << ", max: " << max(D) << endl;
}


void leastSquaresVolPointInterpolation::
makeMirrorPlaneTransformation() const
{
    if (debug)
    {
        Info<< "leastSquaresVolPointInterpolation::"
            << "makeMirrorPlaneTransformation() : "
            << "constructing mirror plane normals and transformation tensors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (mirrorPlaneTransformationPtr_)
    {
        FatalErrorIn
        (
            "leastSquaresVolPointInterpolation::"
            "makeMirrorPlaneTransformation() const"
        )
            << "Mirror plane normals and transformation tensors already exist"
                << abort(FatalError);
    }

    mirrorPlaneTransformationPtr_ =
        new List<Tuple2<vector, tensor> >
        (
            mesh().points().size(),
            Tuple2<vector, tensor>(vector::zero, tensor::zero)
        );
    List<Tuple2<vector, tensor> >& mirrorPlaneTransformation =
        *mirrorPlaneTransformationPtr_;

//     mirrorPlaneTransformationPtr_ = new Map<Tuple2<vector, tensor> >();
//     Map<Tuple2<vector, tensor> >& mirrorPlaneTransformation =
//         *mirrorPlaneTransformationPtr_;


    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            (
                mesh().boundaryMesh()[patchI].type()
             == emptyPolyPatch::typeName
            )
        )
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vectorField& pointNormals =
                mesh().boundaryMesh()[patchI].pointNormals();

            forAll(meshPoints, pointI)
            {
                mirrorPlaneTransformation[meshPoints[pointI]] =
                    Tuple2<vector, tensor>
                    (
                        pointNormals[pointI],
                        I
                    );

//                 mirrorPlaneTransformation.insert
//                 (
//                     meshPoints[pointI],
//                     Tuple2<vector, tensor>
//                     (
//                         pointNormals[pointI],
//                         I
//                     )
//                 );
            }
        }
        else if
        (
            mesh().boundaryMesh()[patchI].type()
         == wedgePolyPatch::typeName
        )
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vectorField& pointNormals =
                mesh().boundaryMesh()[patchI].pointNormals();

            const wedgePolyPatch& wedge =
                refCast<const wedgePolyPatch>(mesh().boundaryMesh()[patchI]);

            forAll(meshPoints, pointI)
            {
                if (!pointAxisEdges().found(meshPoints[pointI]))
                {
                    mirrorPlaneTransformation[meshPoints[pointI]] =
                        Tuple2<vector, tensor>
                        (
                            pointNormals[pointI],
                            wedge.cellT()
                        );
                }
//                 mirrorPlaneTransformation.insert
//                 (
//                     meshPoints[pointI],
//                     Tuple2<vector, tensor>
//                     (
//                         pointNormals[pointI],
//                         wedge.cellT()
//                     )
//                 );
            }
        }
    }
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

leastSquaresVolPointInterpolation::leastSquaresVolPointInterpolation
(
    const fvMesh& vm
)
:
    MeshObject<fvMesh, leastSquaresVolPointInterpolation>(vm),
    pointBndFacesPtr_(NULL),
    pointCyclicFacesPtr_(NULL),
    pointCyclicBndFacesPtr_(NULL),
    pointCyclicGgiFacesPtr_(NULL),
    pointProcFacesPtr_(NULL),
    axisEdgesPtr_(NULL),
    pointAxisEdgesPtr_(NULL),
//     pointNgbProcBndFaceCentresPtr_(NULL),
    globalPointNgbProcBndFaceCentresPtr_(NULL),
    globalPointNgbProcCellCentresPtr_(NULL),
    procBndFacesPtr_(NULL),
    procBndFaceCentresPtr_(NULL),
    pointProcBndFacesPtr_(NULL),
    procCellsPtr_(NULL),
    pointProcCellsPtr_(NULL),
    procCellCentresPtr_(NULL),
    weightsPtr_(NULL),
    originsPtr_(NULL),
    mirrorPlaneTransformationPtr_(NULL),
    invLsMatrices_(0),
    refLPtr_(NULL)
{
//     Pout << mesh().globalData().sharedPointLabels() << endl;

//     if (Pstream::myProcNo() == 3)
//     {
//         Pout << mesh().globalData().sharedPointLabels() << endl;
//         Pout << mesh().globalData().sharedPointAddr() << endl;
//     }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

leastSquaresVolPointInterpolation::~leastSquaresVolPointInterpolation()
{
    deleteDemandDrivenData(pointBndFacesPtr_);
    deleteDemandDrivenData(pointCyclicFacesPtr_);
    deleteDemandDrivenData(pointCyclicBndFacesPtr_);
    deleteDemandDrivenData(pointCyclicGgiFacesPtr_);
    deleteDemandDrivenData(pointProcFacesPtr_);
    deleteDemandDrivenData(axisEdgesPtr_);
    deleteDemandDrivenData(pointAxisEdgesPtr_);
    deleteDemandDrivenData(procBndFacesPtr_);
    deleteDemandDrivenData(pointProcBndFacesPtr_);
    deleteDemandDrivenData(procCellsPtr_);
    deleteDemandDrivenData(pointProcCellsPtr_);
//     deleteDemandDrivenData(pointNgbProcBndFaceCentresPtr_);

    deleteDemandDrivenData(globalPointNgbProcBndFaceCentresPtr_);
    deleteDemandDrivenData(globalPointNgbProcCellCentresPtr_);
    deleteDemandDrivenData(procBndFaceCentresPtr_);
    deleteDemandDrivenData(procCellCentresPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(originsPtr_);
    deleteDemandDrivenData(mirrorPlaneTransformationPtr_);
//     invLsMatrices_.clear();
    deleteDemandDrivenData(refLPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool leastSquaresVolPointInterpolation::movePoints() const
{
    // Original movePoints function
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(originsPtr_);
    deleteDemandDrivenData(refLPtr_);
    deleteDemandDrivenData(mirrorPlaneTransformationPtr_);
    deleteDemandDrivenData(globalPointNgbProcBndFaceCentresPtr_);
    deleteDemandDrivenData(globalPointNgbProcCellCentresPtr_);
    deleteDemandDrivenData(procCellCentresPtr_);
    deleteDemandDrivenData(procBndFaceCentresPtr_);
    invLsMatrices_.clear();

    return true;
}


bool leastSquaresVolPointInterpolation::updateMesh(const mapPolyMesh&) const
{
    deleteDemandDrivenData(pointBndFacesPtr_);
    deleteDemandDrivenData(pointCyclicFacesPtr_);
    deleteDemandDrivenData(pointCyclicBndFacesPtr_);
    deleteDemandDrivenData(pointCyclicGgiFacesPtr_);
    deleteDemandDrivenData(pointProcFacesPtr_);
    deleteDemandDrivenData(axisEdgesPtr_);
    deleteDemandDrivenData(pointAxisEdgesPtr_);
    deleteDemandDrivenData(globalPointNgbProcBndFaceCentresPtr_);
    deleteDemandDrivenData(globalPointNgbProcCellCentresPtr_);

    deleteDemandDrivenData(procBndFacesPtr_);
    deleteDemandDrivenData(procBndFaceCentresPtr_);
    deleteDemandDrivenData(pointProcBndFacesPtr_);

    deleteDemandDrivenData(procCellsPtr_);
    deleteDemandDrivenData(pointProcCellsPtr_);
    deleteDemandDrivenData(procCellCentresPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(originsPtr_);
    deleteDemandDrivenData(refLPtr_);
    deleteDemandDrivenData(mirrorPlaneTransformationPtr_);
    invLsMatrices_.clear();
    deleteDemandDrivenData(refLPtr_);

    return true;
}

const labelListList& leastSquaresVolPointInterpolation::pointBndFaces() const
{
    if (!pointBndFacesPtr_)
    {
        makePointFaces();
    }

    return *pointBndFacesPtr_;
}

const labelListList& leastSquaresVolPointInterpolation
::pointCyclicFaces() const
{
    if (!pointCyclicFacesPtr_)
    {
        makePointFaces();
    }

    return *pointCyclicFacesPtr_;
}

const labelListList& leastSquaresVolPointInterpolation
::pointCyclicBndFaces() const
{
    if (!pointCyclicBndFacesPtr_)
    {
        makePointFaces();
    }

    return *pointCyclicBndFacesPtr_;
}

const labelListList& leastSquaresVolPointInterpolation
::pointCyclicGgiFaces() const
{
    if (!pointCyclicGgiFacesPtr_)
    {
        makePointFaces();
    }

    return *pointCyclicGgiFacesPtr_;
}

const labelListList& leastSquaresVolPointInterpolation
::pointCyclicGgiBndFaces() const
{
    if (!pointCyclicGgiBndFacesPtr_)
    {
        makePointFaces();
    }

    return *pointCyclicGgiBndFacesPtr_;
}

const labelList& leastSquaresVolPointInterpolation
::axisEdges() const
{
    if (!axisEdgesPtr_)
    {
        makeAxisEdges();
    }

    return *axisEdgesPtr_;
}

const Map<labelList>&
leastSquaresVolPointInterpolation::pointAxisEdges() const
{
    if (!pointAxisEdgesPtr_)
    {
        makePointAxisEdges();
    }

    return *pointAxisEdgesPtr_;
}

// const Map<Field<vector> >&
// leastSquaresVolPointInterpolation::pointNgbProcBndFaceCentres() const
// {
//     if (!pointNgbProcBndFaceCentresPtr_)
//     {
//         makePointNgbProcBndFaceCentres();
//     }

//     return *pointNgbProcBndFaceCentresPtr_;
// }

const Map<Field<vector> >&
leastSquaresVolPointInterpolation::globalPointNgbProcBndFaceCentres() const
{
    if (!globalPointNgbProcBndFaceCentresPtr_)
    {
        makeGlobalPointNgbProcBndFaceCentres();
    }

    return *globalPointNgbProcBndFaceCentresPtr_;
}

const Map<Field<vector> >&
leastSquaresVolPointInterpolation::globalPointNgbProcCellCentres() const
{
    if (!globalPointNgbProcCellCentresPtr_)
    {
        makeGlobalPointNgbProcCellCentres();
    }

    return *globalPointNgbProcCellCentresPtr_;
}

const labelListList& leastSquaresVolPointInterpolation::pointProcFaces() const
{
    if (!pointProcFacesPtr_)
    {
        makePointFaces();
    }

    return *pointProcFacesPtr_;
}

const labelListList& leastSquaresVolPointInterpolation::procBndFaces() const
{
    if (!procBndFacesPtr_)
    {
        makeProcBndFaces();
    }

    return *procBndFacesPtr_;
}

const FieldField<Field, vector>&
leastSquaresVolPointInterpolation::procBndFaceCentres() const
{
    if (!procBndFaceCentresPtr_)
    {
        makeProcBndFaceCentres();
    }

    return *procBndFaceCentresPtr_;
}

const List<List<labelPair> >& leastSquaresVolPointInterpolation::
pointProcBndFaces() const
{
    if (!pointProcBndFacesPtr_)
    {
        makeProcBndFaces();
    }

    return *pointProcBndFacesPtr_;
}

const labelListList& leastSquaresVolPointInterpolation::procCells() const
{
    if (!procCellsPtr_)
    {
        makeProcCells();
    }

    return *procCellsPtr_;
}

const Map<List<labelPair> >& leastSquaresVolPointInterpolation::
pointProcCells() const
{
    if (!pointProcCellsPtr_)
    {
        makeProcCells();
    }

    return *pointProcCellsPtr_;
}

const FieldField<Field, vector>&
leastSquaresVolPointInterpolation::procCellCentres() const
{
    if (!procCellCentresPtr_)
    {
        makeProcCellCentres();
    }

    return *procCellCentresPtr_;
}

const FieldField<Field, scalar>&
leastSquaresVolPointInterpolation::weights() const
{
    if (!weightsPtr_)
    {
        makeWeights();
    }

    return *weightsPtr_;
}

const vectorField& leastSquaresVolPointInterpolation::origins() const
{
    if (!originsPtr_)
    {
        makeOrigins();
    }

    return *originsPtr_;
}

const scalarField& leastSquaresVolPointInterpolation::refL() const
{
    if (!refLPtr_)
    {
        makeOrigins();
    }

    return *refLPtr_;
}

const List<Tuple2<vector, tensor> >& leastSquaresVolPointInterpolation::
mirrorPlaneTransformation() const
{
    if (!mirrorPlaneTransformationPtr_)
    {
        makeMirrorPlaneTransformation();
    }

    return *mirrorPlaneTransformationPtr_;
}

const PtrList<scalarRectangularMatrix>&
leastSquaresVolPointInterpolation::invLsMatrices() const
{
    label size = invLsMatrices_.size();

    reduce(size, maxOp<label>());

    if (size == 0)
    {
        makeInvLsMatrices();
    }

    return invLsMatrices_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
