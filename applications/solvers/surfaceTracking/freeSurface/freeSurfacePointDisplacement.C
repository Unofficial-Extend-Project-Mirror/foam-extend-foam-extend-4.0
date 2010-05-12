/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Z. Tukovic and H. Jasak
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

    const edgeList& edges = aMesh().patch().edges();

    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        tensor M = tensor::zero;

        vector S = vector::zero;


        scalarField w(curPointFaces.size(), 0.0);

        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            w[faceI] = 1.0/mag
            (
                controlPoints()[curFace]
              - points[curPoint]
            );
        }

        w /= sum(w);


        forAll (curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            M = M + sqr(w[faceI])*sqr(controlPoints()[curFace]);

            S += sqr(w[faceI])*controlPoints()[curFace];
        }

        vector N = inv(M) & S;

        N /= mag(N);

        scalar p = (S&N)/sum(sqr(w));

        displacement[curPoint] = 
            pointsDisplacementDir()[curPoint]*
            (p - (points[curPoint]&N))/
            (pointsDisplacementDir()[curPoint]&N);
    }



    // Calculate displacement of points which belonge to empty patches

    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    vectorField pointNormals = aMesh().pointAreaNormals();
            
    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];
        
        if (motionPointsMask()[curPoint])
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
                    vector nE = 
                        pointNormals[edges[curEdge].start()]
                      + pointNormals[edges[curEdge].end()];

                    nE /= mag(nE);

                    vector eP =
                        controlPoints()[edgeFaces[curEdge][0]]
                      - edges[curEdge].centre(points);

                    mirrorPoints[++counter] =
                        edges[curEdge].centre(points)
                      + ((2.0*nE*nE - I)&eP);
                }
            }


            // Calculating LS plane interpolation
            const labelList& curPointFaces = pointFaces[curPoint];
        
            tensor M = tensor::zero;
                
            vector S = vector::zero;
                
            scalarField w(curPointFaces.size() + 2, 0.0);
                
            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];
                    
                w[faceI] = 1.0/mag
                    (
                        controlPoints()[curFace]
                      - points[curPoint]
                    );
            }

            forAll (mirrorPoints, pI)
            {
                w[curPointFaces.size() + pI] = 1.0/mag
                    (
                        mirrorPoints[pI]
                      - points[curPoint]
                    );
            }

            w /= sum(w);


            forAll (curPointFaces, faceI)
            {
                label curFace = curPointFaces[faceI];

                M = M + sqr(w[faceI])*sqr(controlPoints()[curFace]);

                S += sqr(w[faceI])*controlPoints()[curFace];
            }

            forAll (mirrorPoints, pI)
            {
                M = M + sqr(w[curPointFaces.size()+pI])*sqr(mirrorPoints[pI]);
            
                S += sqr(w[curPointFaces.size()+pI])*mirrorPoints[pI];
            }


            vector N = inv(M)&S;

            N /= mag(N);

            scalar p = (S&N)/sum(sqr(w));

            displacement[curPoint] = 
                pointsDisplacementDir()[curPoint]*
                (p - (points[curPoint]&N))/
                (pointsDisplacementDir()[curPoint]&N);
        }
    }
            

    forAll(aMesh().boundary(), patchI)
    {
        bool fixedPatch = false;

        forAll(fixedFreeSurfacePatches_, fpI)
        {
            label fixedPatchID = aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[fpI]
            );

            if (fixedPatchID == patchI)
            {
                fixedPatch = true;
            }
        }

        if 
        (
            ( 
                aMesh().boundary()[patchI].type()
             != emptyFaPatch::typeName 
            )
            && !fixedPatch
        )
        {
            labelList patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            labelListList patchPointEdges = 
                aMesh().boundary()[patchI].pointEdges();

            unallocLabelList patchEdgeFaces = 
                aMesh().boundary()[patchI].edgeFaces();

            forAll(patchPoints, pointI)
            {
                forAll(patchPointEdges[pointI], edgeI)
                {
                    label curEdge = patchPointEdges[pointI][edgeI];

                    displacement[patchPoints[pointI]] += 
                        pointsDisplacementDir()[patchPoints[pointI]]*
                        deltaH[patchEdgeFaces[curEdge]];
                }

                displacement[patchPoints[pointI]] /= 
                    patchPointEdges[pointI].size();
            }
        }
    }


    return tdisplacement;
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
