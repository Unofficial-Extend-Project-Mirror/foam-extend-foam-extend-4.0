/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/


#include "freeSurface.H"
#include "primitivePatchInterpolation.H"
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "PstreamCombineReduceOps.H"
#include "coordinateSystem.H"
#include "scalarMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


tmp<vectorField> freeSurface::pointDisplacement(const scalarField& deltaH)
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    controlPoints() += facesDisplacementDir()*deltaH;

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );
    
    vectorField& displacement = tdisplacement();


    // Calculate displacement of internal points
    const vectorField& pointNormals = aMesh().pointAreaNormals();
    const edgeList& edges = aMesh().patch().edges();
    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        vectorField lsPoints(curPointFaces.size(), vector::zero);

        for (label i=0; i<curPointFaces.size(); i++)
        {
            label curFace = curPointFaces[i];
            
            lsPoints[i] = controlPoints()[curFace];
        }
       
        vectorField pointAndNormal = 
            lsPlanePointAndNormal
            (
                lsPoints, 
                points[curPoint], 
                pointNormals[curPoint]
            );

        vector& P = pointAndNormal[0];
        vector& N = pointAndNormal[1];

        displacement[curPoint] = 
            pointsDisplacementDir()[curPoint]
           *((P - points[curPoint])&N)
           /(pointsDisplacementDir()[curPoint]&N);
    }


    // Mirror control points
    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

    forAll(patchMirrorPoints, patchI)
    {
        patchMirrorPoints.set
        (
            patchI,
            new vectorField
            (
                aMesh().boundary()[patchI].faPatch::size(), 
                vector::zero
            )
        );

        vectorField N = 
            aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();
        
        const labelList peFaces = 
            labelList::subList
            (
                aMesh().edgeOwner(),
                aMesh().boundary()[patchI].faPatch::size(),
                aMesh().boundary()[patchI].start()
            );

        const labelList& pEdges = aMesh().boundary()[patchI];

        vectorField peCentres(pEdges.size(), vector::zero);
        forAll(peCentres, edgeI)
        {
            peCentres[edgeI] = 
                edges[pEdges[edgeI]].centre(points);
        }

        vectorField delta =
            vectorField(controlPoints(), peFaces)
          - peCentres;

        patchMirrorPoints[patchI] = 
            peCentres + ((I - 2*N*N)&delta);
    }


    // Calculate displacement of boundary points 
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];
        
        if (motionPointsMask()[curPoint] == 1)
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;
        
            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if(edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;
                    label edgeID = -1;
                    forAll(aMesh().boundary(), patchI)
                    {
                        const labelList& pEdges =
                            aMesh().boundary()[patchI];
                        label index = findIndex(pEdges, curEdge);
                        if (index != -1)
                        {
                            patchID = patchI;
                            edgeID = index;
                            break;
                        }
                    }

                    mirrorPoints[++counter] = 
                        patchMirrorPoints[patchID][edgeID];
                }
            }

            // Calculating LS plane fit
            const labelList& curPointFaces = pointFaces[curPoint];
        
            vectorField lsPoints
            (
                curPointFaces.size() + mirrorPoints.size(), 
                vector::zero
            );
                
            counter = -1;

            for (label i=0; i<curPointFaces.size(); i++)
            {
                label curFace = curPointFaces[i];

                lsPoints[++counter] = controlPoints()[curFace];
            }

            for (label i=0; i<mirrorPoints.size(); i++)
            {
                lsPoints[++counter] = mirrorPoints[i];
            }

            vectorField pointAndNormal = 
                lsPlanePointAndNormal
                (
                    lsPoints, 
                    points[curPoint], 
                    pointNormals[curPoint]
                );

            vector& P = pointAndNormal[0];
            vector& N = pointAndNormal[1];

            displacement[curPoint] = 
                pointsDisplacementDir()[curPoint]
               *((P - points[curPoint])&N)
               /(pointsDisplacementDir()[curPoint]&N);
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type() 
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                label axisPoint = wedgePatch.axisPoint();
                
                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                       &(
                            controlPoints()[pointFaces[axisPoint][0]]
                          - points[axisPoint]
                        )
                    );
            }
        }
    }


    // Calculate displacement of processor patch points
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type() 
         == processorFaPatch::typeName
        )
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPointLabels =
                procPatch.pointLabels();

            FieldField<Field, vector> lsPoints(patchPointLabels.size());
            forAll(lsPoints, pointI)
            {
                lsPoints.set(pointI, new vectorField(0, vector::zero));
            }

            const labelList& nonGlobalPatchPoints = 
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint = 
                    patchPointLabels[curPatchPoint];

                const labelList& curPointFaces = pointFaces[curPoint];

                lsPoints[curPatchPoint].setSize(curPointFaces.size());

                forAll(curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
                }

#               include "boundaryProcessorFaPatchPoints.H"
            }

            scalar lsPointsSize = 0;
            forAll(lsPoints, pointI)
            {
                lsPointsSize +=
                    2*lsPoints[pointI].size()*sizeof(vector);
            }
            
            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(), 
                    lsPointsSize
                );

                toNeighbProc << lsPoints;
            }

            FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    lsPointsSize
                );

                fromNeighbProc >> ngbLsPoints;
            }

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

                vectorField allLsPoints
                (
                    lsPoints[curPatchPoint].size()
                  + ngbLsPoints[curNgbPoint].size(),
                    vector::zero
                );

                label counter = -1;
                forAll(lsPoints[curPatchPoint], pointI)
                {
                    allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
                }
                forAll(ngbLsPoints[curNgbPoint], pointI)
                {
                    allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
                }

                vectorField pointAndNormal =
                    lsPlanePointAndNormal
                    (
                        allLsPoints,
                        points[curPoint],
                        pointNormals[curPoint]
                    );

                vector& P = pointAndNormal[0];
                vector& N = pointAndNormal[1];

                if (motionPointsMask()[curPoint] != 0)
                {
                    displacement[curPoint] =
                        pointsDisplacementDir()[curPoint]
                       *((P - points[curPoint])&N)
                       /(pointsDisplacementDir()[curPoint]&N);
                }
            }
        }
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        const labelList& addr = aMesh().globalData().sharedPointAddr();

        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
        {
            List<List<vector> > procLsPoints(Pstream::nProcs());

            label curSharedPointIndex = findIndex(addr, k);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelList& curPointFaces = pointFaces[curPoint];

                procLsPoints[Pstream::myProcNo()] = 
                    List<vector>(curPointFaces.size());

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    procLsPoints[Pstream::myProcNo()][faceI] = 
                        controlPoints()[curFace];
                }
            }

            Pstream::gatherList(procLsPoints);
            Pstream::scatterList(procLsPoints);
                
            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];
                
                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, vector::zero);

                label counter = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        allPoints[counter++] =
                            procLsPoints[procI][pointI];
                    }
                }

                vectorField pointAndNormal = 
                    lsPlanePointAndNormal
                    (
                        allPoints, 
                        points[curPoint], 
                        pointNormals[curPoint]
                    );

                const vector& P = pointAndNormal[0];
                const vector& N = pointAndNormal[1];

                displacement[curPoint] = 
                    pointsDisplacementDir()[curPoint]
                   *((P - points[curPoint])&N)
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }

    return tdisplacement;
}


tmp<vectorField> freeSurface::lsPlanePointAndNormal
(
    const vectorField& points,
    const vector& origin,
    const vector& axis
) const
{
    // LS in local CS
    vector dir = (points[0] - origin);
    dir -= axis*(axis&dir);
    dir /= mag(dir);
    coordinateSystem cs("cs", origin, axis, dir);

    vectorField localPoints = cs.localPosition(points);
    scalarField W = 1.0/(mag(points - origin) + SMALL);

    scalarRectangularMatrix M
    (
        points.size(),
        3,
        0.0
    );

    for (label i=0; i<localPoints.size(); i++)
    {
        M[i][0] = localPoints[i].x();
        M[i][1] = localPoints[i].y();
        M[i][2] = 1.0;
    }

    scalarSquareMatrix MtM(3, 0.0);
    for (label i = 0; i < MtM.n(); i++)
    {
        for (label j = 0; j < MtM.m(); j++)
        {
            for (label k = 0; k < M.n(); k++)
            {
                MtM[i][j] += M[k][i]*M[k][j]*W[k];
            }
        }
    }

    scalarField MtR(3, 0);
    for (label i = 0; i < MtR.size(); i++)
    {
        for (label j = 0; j < M.n(); j++)
        {
            MtR[i] += M[j][i]*localPoints[j].z()*W[j];
        }
    }

    scalarSquareMatrix::LUsolve(MtM, MtR);

    vector n0 = vector(-MtR[0], -MtR[1], 1);
    n0 = cs.globalVector(n0);
    n0 /= mag(n0);

    vector p0 = vector(0, 0, MtR[2]);
    p0 = cs.globalPosition(p0);

    tmp<vectorField> pointAndNormal
    (
        new vectorField(2, vector::zero)
    );

    pointAndNormal()[0] = p0;
    pointAndNormal()[1] = n0;

    return pointAndNormal;
}


// tmp<vectorField> freeSurface::pointDisplacementSM()
// {
//     const pointField& points = aMesh().patch().localPoints();
//     const labelListList& pointFaces = aMesh().patch().pointFaces();


//     tmp<vectorField> tdisplacement
//     (
//         new vectorField
//         (
//             points.size(),
//             vector::zero
//         )
//     );
    
//     vectorField& displacement = tdisplacement();


//     forAll (pointFaces, pointI)
//     {
//         scalar weightsSum = 0.0;

//         const labelList& curPointFaces = pointFaces[pointI];

//         forAll (curPointFaces, faceI)
//         {
//             label curFace = curPointFaces[faceI];

//             scalar weight = 1.0/mag
//             (
//                 points[pointI]
//               - controlPoints()[curFace]
//             );

//             displacement[pointI] += weight*controlPoints()[curFace];
            
//             weightsSum += weight;
//         }

//         displacement[pointI] /= weightsSum;
        
//         displacement[pointI] -= points[pointI];
//     }


//     displacement = motionPointsMask()*
//         (pointsDisplacementDir()&displacement)*
//         pointsDisplacementDir();


//     return tdisplacement;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
