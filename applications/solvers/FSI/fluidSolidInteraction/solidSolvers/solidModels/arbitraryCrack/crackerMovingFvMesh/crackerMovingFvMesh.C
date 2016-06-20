/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "crackerMovingFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "tetMotionSolver.H"
// #include "laplaceTetMotionSolver.H"
#include "fixedValueTetPolyPatchFields.H"
#include "transformField.H"

#include "eig3.H"

#include "twoDPointCorrector.H"
#include "emptyPolyPatch.H"
#include "unsTotalLagrangianSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(crackerMovingFvMesh, 0);
    addToRunTimeSelectionTable(topoChangerFvMesh, crackerMovingFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::crackerMovingFvMesh::crackerMovingFvMesh(const IOobject& io)
:
    crackerFvMesh(typeName, io),
    motionPtr_
    (
        new laplaceTetMotionSolver
        (
            *this,
            coeffDict().lookup("notchPatchName")
        )
    ),
    notchPatchName_
    (
        coeffDict().lookup("notchPatchName")
    ),
    notchPatchIndex_(-1),
    crackTipEdgeIndex_(-1),
    maxCosFaceIndex_(-1),
    springMeshMotion_(coeffDict().lookup("springMeshMotion")),
    deltaT_(readScalar(coeffDict().lookup("deltaT"))),
    dpTol_(readScalar(coeffDict().lookup("dpTol"))),
    maxIters_(readInt(coeffDict().lookup("maxIters"))),
    fixedPoints_(),
    oldPoints_(points())
{
    notchPatchIndex_ = boundaryMesh().findPatchID(notchPatchName_);

    if(notchPatchIndex_ < 0)
    {
        FatalErrorIn
        (
            "crackerMovingFvMesh::crackerMovingFvMesh(const IOobject& io)"
        )
            << "Can't find patch: " << notchPatchName_
                << exit(FatalError);
    }

//     Info << coeffDict().lookup("springMeshMotion") << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::crackerMovingFvMesh::~crackerMovingFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::crackerMovingFvMesh::updateCrackTip()
{
    // Find crack tip edge
    crackTipEdgeIndex_ = -1;
    vector crackTipEdgeNormal_ = vector::zero;
    {
        label patchIndex = notchPatchIndex_;

        if (boundaryMesh()[crackPatchIndex()].size())
        {
            Info << "crack patch: " <<  crackPatchIndex() << ", "
                << boundaryMesh()[crackPatchIndex()].size() << endl;
            patchIndex = crackPatchIndex();
        }

        const labelList& mEdges =
            boundaryMesh()[patchIndex].meshEdges();

        const labelListList& eFaces =
            boundaryMesh()[patchIndex].edgeFaces();

//         const labelListList& fFaces =
//             boundaryMesh()[patchIndex].faceFaces();

//         const labelListList& pFaces =
//             boundaryMesh()[patchIndex].pointFaces();

        const Field<vector>& fNormals =
            boundaryMesh()[patchIndex].faceNormals();

//         const vectorField::subField fCentres =
        const Field<vector>& fCentres =
            boundaryMesh()[patchIndex].faceCentres();

//         Info << fCentres << endl;
//         Info << fNormals << endl;
//         Info << boundaryMesh()[patchIndex].edges() << endl;
//         Info << eFaces << endl;
//         Info << fFaces << endl;

        forAll(mEdges, edgeI)
        {
            if (eFaces[edgeI].size() == 2)
            {
                if
                (
                    (
                        fNormals[eFaces[edgeI][0]]
                      & fNormals[eFaces[edgeI][1]]
                    ) < -SMALL
                )
                {
                    crackTipEdgeIndex_ = mEdges[edgeI];

                    vector eTilda =
                        edges()[crackTipEdgeIndex_].vec(points());
                    eTilda /= mag(eTilda);

                    crackTipEdgeNormal_ =
                        edges()[crackTipEdgeIndex_].centre(points())
                      - fCentres[eFaces[edgeI][0]];
                    crackTipEdgeNormal_ -=
                        eTilda*(eTilda & crackTipEdgeNormal_);
                    crackTipEdgeNormal_ /= mag(crackTipEdgeNormal_) + SMALL;

//                     Info << fNormals[eFaces[edgeI][0]] << ", "
//                         << fNormals[eFaces[edgeI][1]] << endl;

//                     Info << fCentres[eFaces[edgeI][0]] << ", "
//                         << fCentres[eFaces[edgeI][1]] << endl;

                    break;
                }
            }
        }

        Info << "Crack tip edge centre: "
            << edges()[crackTipEdgeIndex_].centre(points()) << endl;
    }


    if (crackTipEdgeIndex_ != -1)
    {
//         const edge& crackTipEdge = edges()[crackTipEdgeIndex_];

        vector crackTipEdgeCentre =
            edges()[crackTipEdgeIndex_].centre(points());

        const labelList& crackTipEdgeFaces =
            edgeFaces()[crackTipEdgeIndex_];

        const symmTensorField& sigma =
            this->lookupObject<surfaceSymmTensorField>("sigmaf")
           .internalField();

        vector crackTipEdgePrincipalDirection = vector::zero;
        {
            // Looking up solid solver
            const solidSolver& stress =
                this->lookupObject<solidSolver>
                (
                    "solidProperties"
                );

            const solidSolvers::unsTotalLagrangianSolid& tlStress =
                refCast<const solidSolvers::unsTotalLagrangianSolid>
                (
                    stress
                );

            pointSymmTensorField pointSigma
            (
                IOobject
                (
                    "pointSigam",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                tlStress.pMesh(),
                dimensioned<symmTensor>
                (
                    "0",
                    tlStress.sigma().dimensions(),
                    symmTensor::zero
                )
            );

            for (direction cmpt = 0; cmpt < symmTensor::nComponents; cmpt++)
            {
                volScalarField cmptSigma = tlStress.sigma().component(cmpt);

                pointScalarField cmptPointSigma
                (
                    IOobject
                    (
                        "cmptPointSigma",
                        time().timeName(),
                        *this,
                        IOobject::NO_READ
                    ),
                    tlStress.pMesh(),
                    dimensioned<scalar>
                    (
                        "0",
                        tlStress.sigma().dimensions(),
                        0
                    )
                );

                tlStress.volToPoint().interpolate(cmptSigma, cmptPointSigma);

                pointSigma.internalField().replace
                (
                    cmpt,
                    cmptPointSigma.internalField()
                );
            }

            const edge& crackTipEdge = edges()[crackTipEdgeIndex_];
            label startPoint = crackTipEdge.start();
            label endPoint = crackTipEdge.end();

            tensor V = tensor::zero;
            vector d = vector::zero;

            eigen_decomposition(pointSigma[startPoint], V, d);
            crackTipEdgePrincipalDirection +=
                vector(V.zx(), V.zy(), V.zz());

            eigen_decomposition(pointSigma[endPoint], V, d);
            crackTipEdgePrincipalDirection +=
                vector(V.zx(), V.zy(), V.zz());

            crackTipEdgePrincipalDirection /=
                mag(crackTipEdgePrincipalDirection);
        }

        vectorField n = faceAreas();
        n /= mag(n);

        scalar maxCos = 0;
        maxCosFaceIndex_ = -1;
        vector maxCosPrincipalDirection = vector::zero;

        forAll(crackTipEdgeFaces, faceI)
        {
            label curFace = crackTipEdgeFaces[faceI];

            if
            (
                isInternalFace(curFace)
             && (
                    crackTipEdgeNormal_
                  & (faceCentres()[curFace] - crackTipEdgeCentre)
                ) > SMALL
             && mag(sigma[curFace]) > SMALL
            )
            {
                tensor V = tensor::zero;
                vector d = vector::zero;
                eigen_decomposition(sigma[curFace], V, d);

//                     scalar principalStress = d.z();
                vector principalDirection(V.zx(), V.zy(), V.zz());

//                 scalar cos = mag(principalDirection & n[curFace]);
                scalar cos = mag(crackTipEdgePrincipalDirection & n[curFace]);
                if (cos > maxCos)
                {
                    maxCos = cos;
                    maxCosFaceIndex_ = curFace;
                    maxCosPrincipalDirection = principalDirection;
                }
            }
        }

        Info << "Crack tip face centre: "
            << faceCentres()[maxCosFaceIndex_] << ", "
            << n[maxCosFaceIndex_] << endl;
    }

    // Store current points
    oldPoints_ = points();

    return true;
}

bool Foam::crackerMovingFvMesh::smoothMesh()
{
    Info << "Smooth mesh" << endl;

    vectorField newPoints = points();

    fixedPoints_.clear();

    smooth(newPoints);

    twoDPointCorrector twoDPointCorrector(*this);
    twoDPointCorrector.correctPoints(newPoints);

    fvMesh::movePoints(newPoints);

    return true;
}

bool Foam::crackerMovingFvMesh::update()
{
//     Info << "springMeshMotion = " << springMeshMotion_ << endl;

    if (crackerFvMesh::update())
    {
        motionPtr_->updateMesh(topoChangeMap());

        return true;
    }
    else // Mesh motion
    {
        if (maxCosFaceIndex_ != -1)
        {
            const symmTensorField& sigma =
                this->lookupObject<surfaceSymmTensorField>("sigmaf")
               .internalField();

            tensor V = tensor::zero;
            vector d = vector::zero;
            eigen_decomposition(sigma[maxCosFaceIndex_], V, d);

            vector maxCosPrincipalDirection(V.zx(), V.zy(), V.zz());

            vector crackTipEdgePrincipalDirection = vector::zero;
            {
                // Looking up solid solver
                const solidSolver& stress =
                    this->lookupObject<solidSolver>
                    (
                        "solidProperties"
                    );

                const solidSolvers::unsTotalLagrangianSolid& tlStress =
                    refCast<const solidSolvers::unsTotalLagrangianSolid>
                    (
                        stress
                    );

                pointSymmTensorField pointSigma
                (
                    IOobject
                    (
                        "pointSigam",
                        time().timeName(),
                        *this,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tlStress.pMesh(),
                    dimensioned<symmTensor>
                    (
                        "0",
                        tlStress.sigma().dimensions(),
                        symmTensor::zero
                    )
                );

                for
                (
                    direction cmpt = 0;
                    cmpt < symmTensor::nComponents;
                    cmpt++
                )
                {
                    volScalarField cmptSigma =
                        tlStress.sigma().component(cmpt);

                    pointScalarField cmptPointSigma
                    (
                        IOobject
                        (
                            "cmptPointSigma",
                            time().timeName(),
                            *this,
                            IOobject::NO_READ
                        ),
                        tlStress.pMesh(),
                        dimensioned<scalar>
                        (
                            "0",
                            tlStress.sigma().dimensions(),
                            0
                        )
                    );

                    tlStress.volToPoint().interpolate
                    (
                        cmptSigma,
                        cmptPointSigma
                    );

                    pointSigma.internalField().replace
                    (
                        cmpt,
                        cmptPointSigma.internalField()
                    );
                }

                const edge& crackTipEdge = edges()[crackTipEdgeIndex_];
                label startPoint = crackTipEdge.start();
                label endPoint = crackTipEdge.end();

                tensor V = tensor::zero;
                vector d = vector::zero;

                eigen_decomposition(pointSigma[startPoint], V, d);
                crackTipEdgePrincipalDirection +=
                    vector(V.zx(), V.zy(), V.zz());

                eigen_decomposition(pointSigma[endPoint], V, d);
                crackTipEdgePrincipalDirection +=
                    vector(V.zx(), V.zy(), V.zz());

                crackTipEdgePrincipalDirection /=
                    mag(crackTipEdgePrincipalDirection);
            }

            maxCosPrincipalDirection = crackTipEdgePrincipalDirection;

            vector n = faceAreas()[maxCosFaceIndex_];
            n /= mag(n);

            Info << "principal vector: "
                << maxCosPrincipalDirection << endl;

            if ( (maxCosPrincipalDirection & n) < -SMALL)
            {
                maxCosPrincipalDirection *= -1;
            }

            vector relaxedNormal = n + maxCosPrincipalDirection;
            relaxedNormal /= mag(relaxedNormal);

            relaxedNormal = n;

            scalar maxCos = mag(maxCosPrincipalDirection & relaxedNormal);

            tensor R =
                rotationTensor
                (
                    relaxedNormal,
                    maxCosPrincipalDirection
                );

            Info << crackTipEdgeIndex_ << ", "
                << maxCosFaceIndex_ << ", "
                << faceCentres()[maxCosFaceIndex_] << ", "
                << n << ", "
                << maxCosPrincipalDirection << ", "
                << maxCos << endl;

            if (mag(mag(maxCos) - 1) > 1e-3)
            {
                if (!springMeshMotion_)
                {
                    Info << "FE mesh motion" << endl;

                    motionPtr_->clearConstraints();

                    label curFaceIndex = maxCosFaceIndex_;

                    const labelList& pointLabels = faces()[curFaceIndex];

                    vectorField curOldPoints(oldPoints_, pointLabels);
                    vectorField curPoints(points(), pointLabels);

                    const edge& crackTipEdge = edges()[crackTipEdgeIndex_];

                    vector curRotationOrigin =
                        points()[crackTipEdge.start()];

                    vectorField curVelocities =
                        transform(R, curPoints - curRotationOrigin)
                      + curRotationOrigin
                      - curOldPoints;
                    curVelocities /= time().deltaT().value();

                    forAll(pointLabels, pointI)
                    {
                        motionPtr_->setConstraint
                        (
                            pointLabels[pointI],
                            curVelocities[pointI]
                        );
                    }

                    vector curFaceVelocity =
                        transform
                        (
                            R,
//                             faceCentres()[curFaceIndex]
                            faces()[curFaceIndex].centre(points())
                          - curRotationOrigin
                        )
                      + curRotationOrigin
                      - faces()[curFaceIndex].centre(oldPoints_);
//                       - faceCentres()[curFaceIndex];
                    curFaceVelocity /= time().deltaT().value();

                    motionPtr_->setConstraint
                    (
                        nPoints() + curFaceIndex,
                        curFaceVelocity
                    );

                    fvMesh::movePoints(oldPoints_);
                    fvMesh::movePoints(motionPtr_->newPoints());
                }
                else
                {
                    Info << "Spring analogy mesh motion" << endl;

                    vectorField newPoints = points();

                    fixedPoints_.clear();

                    const labelList& pointLabels = faces()[maxCosFaceIndex_];

                    vectorField curPoints(newPoints, pointLabels);

                    const edge& crackTipEdge = edges()[crackTipEdgeIndex_];

                    vector curRotationOrigin =
                        newPoints[crackTipEdge.start()];

                    vectorField curDisplacements =
                        transform(R, curPoints - curRotationOrigin)
                      + curRotationOrigin
                      - curPoints;

                    forAll(pointLabels, pointI)
                    {
                        fixedPoints_.append(pointLabels[pointI]);
                        newPoints[pointLabels[pointI]] +=
                            curDisplacements[pointI];
                    }

                    smooth(newPoints);

                    twoDPointCorrector twoDPointCorrector(*this);
                    twoDPointCorrector.correctPoints(newPoints);

                    fvMesh::movePoints(newPoints);
                }
            }
        }
    }

    // Mesh motion only - return false
    return false;
}

Foam::scalar Foam::crackerMovingFvMesh::sizeFunction
(
    const vector& p
) const
{
    scalar sf = 1.0;

    return sf;
}

void Foam::crackerMovingFvMesh::smooth(vectorField& points)
{
    Info << "Spring analogy mesh motion" << endl;

    scalar maxRelDisp = 1;
    label iCorr = 0;

    const edgeList& e = edges();

    scalar Fscale_ = 1.2;

    scalar h0_ = 0;
    label n = 0;
    forAll(e, edgeI)
    {
        scalar eL = e[edgeI].mag(points);
        vector eTilda = e[edgeI].vec(points);
        eTilda /= mag(eTilda);

        if (mag(eTilda.z()) < 0.01)
        {
            h0_ += eL;
            n++;
        }
    }
    h0_ /= n;
    Info << "h0 = " << h0_ << endl;

    do
    {
        vectorField Ftot(points.size(), vector::zero);

        vectorField eVec(e.size(), vector::zero);
        scalarField eL(e.size(), 0.0);
        scalarField eSizeFunc(e.size(), 0.0);

        scalar sumL2 = 0;
        scalar sumSizeFunc2 = 0;

        forAll(e, edgeI)
        {
            eVec[edgeI] = e[edgeI].vec(points);
            eL[edgeI] = e[edgeI].mag(points);

            vector eMidPoint = e[edgeI].centre(points);
            eSizeFunc[edgeI] = sizeFunction(eMidPoint);

            vector eTilda = e[edgeI].vec(points);
            eTilda /= mag(eTilda);

            if (mag(eTilda.z()) < 0.01)
            {
                sumL2 += sqr(eL[edgeI]);
                sumSizeFunc2 += sqr(eSizeFunc[edgeI]);
            }
        }

        scalarField eL0 =
            eSizeFunc*Fscale_*sqrt(sumL2/sumSizeFunc2);

        forAll(e, edgeI)
        {
            scalar f = eL0[edgeI] - eL[edgeI];
            if (f<0)
            {
                f = 0;
            }

            vector F = f*eVec[edgeI]/eL[edgeI];

            Ftot[e[edgeI].end()] += F;
            Ftot[e[edgeI].start()] -= F;
        }

        forAll(boundaryMesh(), patchI)
        {
            if
            (
                boundaryMesh()[patchI].type()
             != emptyPolyPatch::typeName
            )
            {
                const labelList& patchPoints =
                    boundaryMesh()[patchI].meshPoints();

                for (label i=0; i<patchPoints.size(); i++)
                {
                    Ftot[patchPoints[i]] = vector::zero;
                }
            }
        }

        forAll (fixedPoints_, i)
        {
            Ftot[fixedPoints_[i]] = vector::zero;
        }

        vectorField delta = deltaT_*Ftot;
        delta.replace(2, 0);

//         Info << delta << endl;

        points += delta;

        maxRelDisp = max(mag(delta))/h0_;
    }
    while(maxRelDisp > dpTol_ && ++iCorr < maxIters_);

    Info << "Num of iterations: " << iCorr
        << ", maxRelDisp: " << maxRelDisp << endl;
}

// bool Foam::crackerMovingFvMesh::update()
// {
//     if (crackerFvMesh::update())
//     {
//         motionPtr_->updateMesh(topoChangeMap());

//         return true;
//     }
//     else // Mesh motion
//     {
//         // Find crack tip edge
//         label crackTipEdgeIndex = -1;
//         vector crackTipEdgeNormal = vector::zero;
//         {
//             // Find initial crack tip

//             label patchIndex = notchPatchIndex_;

//             if (boundaryMesh()[crackPatchIndex()].size())
//             {
//                 patchIndex = crackPatchIndex();
//             }

//             const labelList& meshEdges =
//                 boundaryMesh()[patchIndex].meshEdges();

//             const labelListList& edgeFaces =
//                 boundaryMesh()[patchIndex].edgeFaces();

//             const Field<vector>& faceNormals =
//                 boundaryMesh()[patchIndex].faceNormals();

//             const Field<vector>& faceCentres =
//                 boundaryMesh()[patchIndex].faceCentres();

//             forAll(meshEdges, edgeI)
//             {
//                 if (edgeFaces[edgeI].size() == 2)
//                 {
//                     if
//                     (
//                         (
//                             faceNormals[edgeFaces[edgeI][0]]
//                           & faceNormals[edgeFaces[edgeI][1]]
//                         ) < -SMALL
//                     )
//                     {
//                         crackTipEdgeIndex = meshEdges[edgeI];

//                         vector eTilda =
//                             edges()[crackTipEdgeIndex].vec(points());
//                         eTilda /= mag(eTilda);

//                         crackTipEdgeNormal =
//                             edges()[crackTipEdgeIndex].centre(points())
//                           - faceCentres[edgeFaces[edgeI][0]];
//                         crackTipEdgeNormal -=
//                             eTilda*(eTilda&crackTipEdgeNormal);
//                         crackTipEdgeNormal /= mag(crackTipEdgeNormal) + SMALL;

//                         break;
//                     }
//                 }
//             }
//         }

//         if (crackTipEdgeIndex != -1)
//         {
//             const edge& crackTipEdge = edges()[crackTipEdgeIndex];

//             vector crackTipEdgeCentre =
//                 edges()[crackTipEdgeIndex].centre(points());

//             const labelList& crackTipEdgeFaces =
//                 edgeFaces()[crackTipEdgeIndex];

//             const symmTensorField& sigma =
//                 this->lookupObject<surfaceSymmTensorField>("sigmaf")
//                .internalField();

//             vectorField n = faceAreas();
//             n /= mag(n);

//             scalar maxCos = 0;
//             label maxCosFaceIndex = -1;
//             vector maxCosPrincipalDirection = vector::zero;

//             forAll(crackTipEdgeFaces, faceI)
//             {
//                 label curFace = crackTipEdgeFaces[faceI];

//                 if
//                 (
//                     isInternalFace(curFace)
//                  && (
//                         crackTipEdgeNormal
//                       & (faceCentres()[curFace] - crackTipEdgeCentre)
//                     ) > SMALL
//                  && mag(sigma[curFace]) > SMALL
//                 )
//                 {
//                     tensor V = tensor::zero;
//                     vector d = vector::zero;
//                     eigen_decomposition(sigma[curFace], V, d);

// //                     scalar principalStress = d.z();
//                     vector principalDirection(V.zx(), V.zy(), V.zz());

//                     scalar cos = mag(principalDirection & n[curFace]);
//                     if (cos > maxCos)
//                     {
//                         maxCos = cos;
//                         maxCosFaceIndex = curFace;
//                         maxCosPrincipalDirection = principalDirection;
//                     }
//                 }
//             }

//             if (maxCosFaceIndex != -1)
//             {
//                 Info << "principal vector: "
//                     << maxCosPrincipalDirection << endl;

//                 if ( (maxCosPrincipalDirection & n[maxCosFaceIndex]) < -SMALL)
//                 {
//                     maxCosPrincipalDirection *= -1;
//                 }

//                 tensor R =
//                     rotationTensor
//                     (
//                         n[maxCosFaceIndex],
//                         maxCosPrincipalDirection
//                     );

//                 motionPtr_->clearConstraints();

//                 forAll(crackTipEdgeFaces, faceI)
//                 {
//                     label curFaceIndex = crackTipEdgeFaces[faceI];

//                     if (curFaceIndex == maxCosFaceIndex)
// //                     if (isInternalFace(curFaceIndex))
//                     {
//                         const labelList& pointLabels = faces()[curFaceIndex];

//                         vectorField curPoints(points(), pointLabels);

//                         vector curRotationOrigin =
//                             points()[crackTipEdge.start()];

//                         vectorField curVelocities =
//                             transform(R, curPoints - curRotationOrigin)
//                           + curRotationOrigin
//                           - curPoints;
//                         curVelocities /= time().deltaT().value();

//                         forAll(pointLabels, pointI)
//                         {
// //                             if (mag(curVelocities[pointI]) > SMALL)
//                             {
//                                 motionPtr_->setConstraint
//                                 (
//                                     pointLabels[pointI],
//                                     curVelocities[pointI]
//                                 );
//                             }
//                         }

//                         vector curFaceVelocity =
//                             transform
//                             (
//                                 R,
//                                 faceCentres()[curFaceIndex]
//                               - curRotationOrigin
//                             )
//                           + curRotationOrigin
//                           - faceCentres()[curFaceIndex];
//                         curFaceVelocity /= time().deltaT().value();

// //                         if (mag(curFaceVelocity) > SMALL)
//                         {
//                             motionPtr_->setConstraint
//                             (
//                                 nPoints() + curFaceIndex,
//                                 curFaceVelocity
//                             );
//                         }
//                     }
//                 }

//                 if (mag(mag(maxCos) - 1) > 1e-3)
//                 {
//                     fvMesh::movePoints(motionPtr_->newPoints());
//                 }

//                 Info << crackTipEdgeIndex << ", "
//                     << maxCosFaceIndex << ", "
//                     << faceCentres()[maxCosFaceIndex] << ", "
//                     << n[maxCosFaceIndex] << ", "
//                     << maxCosPrincipalDirection << ", "
//                     << maxCos << endl;

// //                 Info << crackTipEdgeIndex << ", "
// //                     << maxCosFaceIndex << ", "
// //                     << faceCentres()[maxCosFaceIndex] << ", "
// //                     << n[maxCosFaceIndex] << ", "
// //                     << maxCosPrincipalDirection << ", "
// //                     << Foam::acos(maxCos)*180/M_PI << endl;
//             }
//         }
//     }

//     // Mesh motion only - return false
//     return false;
// }

// ************************************************************************* //
