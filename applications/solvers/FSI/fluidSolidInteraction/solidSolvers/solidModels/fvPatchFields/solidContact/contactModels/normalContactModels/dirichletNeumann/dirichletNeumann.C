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

Class
    dirichletNeumann

\*---------------------------------------------------------------------------*/

#include "dirichletNeumann.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "tractionBoundaryGradient.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(dirichletNeumann, 0);
  addToRunTimeSelectionTable(normalContactModel, dirichletNeumann, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void dirichletNeumann::calculateSlaveFaceNormals
// (
//     vectorField& slaveFaceNormals,
//     const nonLinearGeometry::nonLinearType nonLinear,
//     const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
//     const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
//     ExtendedGGIInterpolation<
//         PrimitivePatch< face, Foam::List, pointField >,
//         PrimitivePatch< face, Foam::List, pointField >
//         >* ggiInterpolatorPtr
// )
// {
//     const fvMesh& mesh = mesh_;
//     const label slavePatchIndex = slavePatchID();

//     // dirichletNeumannCalculateSlaveFaceNormals.H
//     if (nonLinear == nonLinearGeometry::OFF)
//     {
//         // undeformed normals
//         // remember that the mesh has not moved, only the global face zone
//         // patches are moved to the deformed position

//         // slaveFaceNormals =
//         //-1*masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
//         //        (
//         //         mesh.boundaryMesh()[masterPatchIndex].faceNormals()
//         //         );
//         // slaveFaceNormals /= mag(slaveFaceNormals);

//         // OPTION 2
//         // Use deformed slave face normals
//         // Get them from global face zone patches
//         // const vectorField& globalSlaveFaceNormals =
//         //slaveFaceZonePatch.faceNormals();
//         // const label slavePatchStart
//         //    = mesh.boundaryMesh()[slavePatchIndex].start();
//         // forAll(slaveFaceNormals, facei)
//         //    {
//         //      slaveFaceNormals[facei] =
//         //        globalSlaveFaceNormals
//         //        [
//         //    mesh.faceZones()[slaveFaceZoneID()].whichFace
//         // (slavePatchStart + facei)
//         //         ];
//         //    }
//         // // make sure they are unity
//         // slaveFaceNormals /= mag(slaveFaceNormals);

//         // OPTION 3
//         // globally interpolate master normals to the slave
//         // then get the local normals

//         // undeformed master normals
//         // const vectorField& actualSlaveFaceNormals =
//         //mesh.boundaryMesh()[slavePatchIndex].faceNormals();
//         // vectorField globalMasterFaceNormals
//         //(masterFaceZonePatch.size(), vector::zero);
//         // const label masterPatchStart
//         //    = mesh.boundaryMesh()[masterPatchIndex].start();
//         // // master undeformed face normals
//         // vectorField masterFaceNormals =
//         //mesh.boundaryMesh()[masterPatchIndex].faceNormals();
//         // // put field into global
//         // forAll(masterFaceNormals, i)
//         //    {
//         //      globalMasterFaceNormals[mesh.faceZones()
//         //[masterFaceZoneID()].whichFace(masterPatchStart + i)] =
//         //        masterFaceNormals[i];
//         //    }

//         // //- exchange parallel data
//         // // sum because each face is only on one proc
//         // reduce(globalMasterFaceNormals, sumOp<vectorField>());

//         // use deformed master normals
//         const vectorField& globalMasterFaceNormals =
//             masterFaceZonePatch.faceNormals();

//         // interpolate master normals to the slave and reverse the direction
//         // we can use inverseDistance or GGI for the interpolation
//         // GGI is much better is we want the actual master normals
//         vectorField globalSlaveFaceNormals
//             (
//                 slaveFaceZonePatch.size(), vector::zero
//             );
//         if (ggiInterpolatorPtr)
//         {
//             //Info << "interpolating master normals with GGI" << endl;
//             globalSlaveFaceNormals =
//                 -ggiInterpolatorPtr->masterToSlave(globalMasterFaceNormals);
//         }
//         else // inverse distance
//         {
//             FatalErrorIn("dirchletNeumannCalculateSlaveFaceNormals.H")
//                 << "patchToPatch disabled, use GGI"
//                     << abort(FatalError);
//             // globalSlaveFaceNormals =
//          //    -masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
//             //     (
//             //         globalMasterFaceNormals
//             //         );
//         }
//         vectorField actualSlaveFaceNormals =
//             slaveFaceZonePatch.faceNormals();

//         // put global back into local
//         const label slavePatchStart
//             = mesh.boundaryMesh()[slavePatchIndex].start();
//         forAll(slaveFaceNormals, facei)
//         {
//             slaveFaceNormals[facei] =
//                 globalSlaveFaceNormals
//                 [
//                     mesh.faceZones()[slaveFaceZoneID()].whichFace
//                     (
//                         slavePatchStart + facei
//                     )
//                 ];

//             if (mag(slaveFaceNormals[facei]) < SMALL)
//             {
//                 // interpolation sometimes does not work for far
//                 // away faces so we will use actual slave normals
//                 // but doesn't really matter as the face is probably
//                 // not in contact
//                 slaveFaceNormals[facei] = actualSlaveFaceNormals[facei];
//             }
//             else
//             {
//                 // make sure they are unity
//                 slaveFaceNormals[facei] /= mag(slaveFaceNormals[facei]);
//             }
//         }
//     }
//     else if
//     (
//         nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
//         || nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
//     )
//     {
//         Info << "nonLinear UL/TL normals" << endl;

//         // a few different ways to calculate deformed normals

//         // OPTION 1
//         // use deformed slave normals calculated using F
//         // calculate deformed normals using deformation gradient
//         // const volTensorField& gradField =
//         // mesh.objectRegistry::lookupObject<volTensorField>
//         //("grad("+fieldName+")");
//         //        // deformation gradient
//         //        tensorField F = gradField.boundaryField()[slavePatchIndex];
//         //        if (fieldName == "DU")
//         //    {
//         //      const volTensorField& gradU =
//      //        mesh.objectRegistry::lookupObject<volTensorField>("grad(U)");
//         //      F += gradU.boundaryField()[slavePatchIndex];
//         //    }
//         //        const tensorField Finv = inv(I + F);
//         // deformed normals
//         //slaveFaceNormals =
//         // J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();


//         // OPTION 2
//         // Use deformed slave face normals
//         // Get them from global face zone patches
//         const vectorField& globalSlaveFaceNormals =
//             slaveFaceZonePatch.faceNormals();
//         const label slavePatchStart
//             = mesh.boundaryMesh()[slavePatchIndex].start();
//         forAll(slaveFaceNormals, facei)
//         {
//             slaveFaceNormals[facei] =
//                 globalSlaveFaceNormals
//                 [
//                     mesh.faceZones()[slaveFaceZoneID()].whichFace
//                     (
//                         slavePatchStart + facei
//                     )
//                 ];
//         }
//         // make sure they are unity
//         slaveFaceNormals /= mag(slaveFaceNormals);


//         // OPTION 3
//         // Use deformed master normals
//         // globally interpolate deformed master normals to the slave
//         // then get the local normals
//         // globally interpolate master normals to the slave
//         // then get the local normals

//         // undeformed master normals
//         //const vectorField& actualSlaveFaceNormals =
//         //mesh.boundaryMesh()[slavePatchIndex].faceNormals();
//         // vectorField globalMasterFaceNormals
//         //(masterFaceZonePatch.size(), vector::zero);
//         // const label masterPatchStart
//         //    = mesh.boundaryMesh()[masterPatchIndex].start();
//         // // master undeformed face normals
//         // vectorField masterFaceNormals =
//         //mesh.boundaryMesh()[masterPatchIndex].faceNormals();
//         // // put field into global
//         // forAll(masterFaceNormals, i)
//         //    {
//         //      globalMasterFaceNormals[mesh.faceZones()
//         //[masterFaceZoneID()].whichFace(masterPatchStart + i)] =
//         //        masterFaceNormals[i];
//         //    }

//         // //- exchange parallel data
//         // // sum because each face is only on one proc
//         // reduce(globalMasterFaceNormals, sumOp<vectorField>());

//         // use deformed master normals
//         // {
//         // const vectorField& globalMasterFaceNormals =
//         //masterFaceZonePatch.faceNormals();

//       // // interpolate master normals to the slave and reverse the direction
//         // // we can use inverseDistance or GGI for the interpolation
//         // // GGI is much better is we want the actual master normals
//         // vectorField globalSlaveFaceNormals
//         //(slaveFaceZonePatch.size(), vector::zero);
//         // if (ggiInterpolatorPtr)
//         //    {
//         //      //Info << "interpolating master normals with GGI" << endl;
//         //      globalSlaveFaceNormals =
//       //        -ggiInterpolatorPtr->masterToSlave(globalMasterFaceNormals);
//         //    }
//         // else // inverse distance
//         //    {
//         //      globalSlaveFaceNormals =
//      //        -masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
//         //        (
//         //         globalMasterFaceNormals
//         //         );
//         //    }
//         // vectorField actualSlaveFaceNormals =
//         //    slaveFaceZonePatch.faceNormals();

//         // // get local normals from global
//         // const label slavePatchStart
//         //    = mesh.boundaryMesh()[slavePatchIndex].start();
//         // forAll(slaveFaceNormals, facei)
//         //    {
//         //      slaveFaceNormals[facei] =
//         //        globalSlaveFaceNormals
//         //        [
//         //         mesh.faceZones()[slaveFaceZoneID()].whichFace
//         //(slavePatchStart + facei)
//         //         ];

//         //      if (mag(slaveFaceNormals[facei]) < SMALL)
//         //        {
//         //          // interpolation sometimes does not work for far
//         //          // away faces so we will use actual slave normals
//         //          // but doesn't really matter as the face is probably
//         //          // not in contact
//         //          slaveFaceNormals[facei] = actualSlaveFaceNormals[facei];
//         //        }
//         //      else
//         //        {
//         //          // make sure they are unity
//         //          slaveFaceNormals[facei] /= mag(slaveFaceNormals[facei]);
//         //        }
//         //    }
//         // }
//     }
//     else if (nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF)
//     {
//         // OPTION 1
//         // use deformed slave normals calculated using F
//         // calculate deformed normals using deformation gradient

//         // Lookup relative deformation gradient
//         // const tensorField& Finv =
//         //     mesh.objectRegistry::lookupObject<volTensorField>
//         //     (
//         //         "relFinv"
//         //     ).boundaryField()[slavePatchIndex];
//         // const scalarField& J =
//         //     mesh.objectRegistry::lookupObject<volScalarField>
//         //     (
//         //         "J"
//         //     ).boundaryField()[slavePatchIndex];

//         // const label masterPatchIndex = masterPatchID();

//         const tensorField& Finv =
//             mesh.objectRegistry::lookupObject<surfaceTensorField>
//             (
//                 "relFinvf"
//             ).boundaryField()[slavePatchIndex];
//         //).boundaryField()[masterPatchIndex];
//         // const scalarField& J =
//         //     mesh.objectRegistry::lookupObject<surfaceScalarField>
//         //     (
//         //         "Jf"
//         //     ).boundaryField()[slavePatchIndex];

//         // Calculate deformed normals by Nanson's formula

//         // if (gMax(mag(slaveFaceNormals_)) < SMALL)
//         // {
//         //     slaveFaceNormals_ =
//       //   //J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
//       //       Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
//         //     slaveFaceNormals_ /= mag(slaveFaceNormals_);
//         // }
//         // const vectorField slaveFaceNormalsPrevIter = slaveFaceNormals_;

//         slaveFaceNormals =
//            //J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
//             Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();

//         slaveFaceNormals /= mag(slaveFaceNormals);
//     }
//     else if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN_KIRCHHOFF)
//     {
//         // OPTION 1
//         // use deformed slave normals calculated using F
//         // calculate deformed normals using deformation gradient

//         // Lookup relative deformation gradient and inverse from the solver
//         const tensorField& Finv =
//             mesh.objectRegistry::lookupObject<volTensorField>
//             (
//                 "Finv"
//             ).boundaryField()[slavePatchIndex];

//         const scalarField& J =
//             mesh.objectRegistry::lookupObject<volScalarField>
//             (
//                 "J"
//             ).boundaryField()[slavePatchIndex];


//         // Calculate deformed normals by Nanson's formula
//         slaveFaceNormals =
//             J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
//     }
//     else if
//     (
//         nonLinear == nonLinearGeometry::DEFORMED_LAGRANGIAN
//     )
//     {
//         //Info << "nonLinear DL normals" << endl;

//         // The main mesh is always at the deformed position for
//         // deformedLagrangian
//      // So we will lookup up the slaveFaceZone from the main mesh as this was
//         // moved using
//         // least squares interpolation which is more accurate than
//         //primitvePatchIntepolation

//         // slave normals
//         // Philipc: fix for parallel - use local patch not global face zone
//         //slaveFaceNormals =
//         //mesh.faceZones()[slaveFaceZoneID()]().faceNormals();
//         slaveFaceNormals = mesh.boundaryMesh()[slavePatchID()].faceNormals();

//         // or interpolate master normals
//         //...
//     }
//     else
//     {
//         FatalError
//             << "dirichletNeumann::correct()" << nl
//                 << "dirichletNeumannCalculateSlaveFaceNormals" << nl
//                 << "nonLinear option " << nonLinear
//                 << " is unknown"
//                 << abort(FatalError);
//     }
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dirichletNeumann::dirichletNeumann
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID,
    const primitiveFacePatch& masterFaceZonePatch,
    const primitiveFacePatch& slaveFaceZonePatch
 )
:
    normalContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZoneID,
        slaveFaceZoneID,
        masterFaceZonePatch,
        slaveFaceZonePatch
    ),
    normalContactModelDict_(dict.subDict(name+"NormalModelDict")),
    mesh_(patch.boundaryMesh().mesh()),
    slaveDisp_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    slavePressure_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    oldSlavePressure_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    zeroSlavePressure_(slavePressure_.size(), vector::zero),
    residual_(1.0),
    maxNorm_(SMALL),
    averagePenetration_(1.0),
    nonLinearType_("off"),
    //touchFraction_(mesh_.boundaryMesh()[slavePatchID].size(), 0.0),
    areaInContact_(mesh_.boundaryMesh()[slavePatchID].size(), 0.0),
    slaveValueFrac_
    (
        mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero
    ),
    oldSlaveValueFrac_
    (
        mesh_.boundaryMesh()[slavePatchID].size(),
        symmTensor::zero
    ),
    limitPenetration_(normalContactModelDict_.lookup("limitPenetration")),
    penetrationLimit_
    (
        readScalar(normalContactModelDict_.lookup("penetrationLimit"))
    ),
    limitPressure_(normalContactModelDict_.lookup("limitPressure")),
    pressureLimit_(readScalar(normalContactModelDict_.lookup("pressureLimit"))),
    correctMissedVertices_
    (
        normalContactModelDict_.lookup("correctMissedVertices")
    ),
    slavePointPointsPtr_(NULL),
    contactGapTol_(readScalar(normalContactModelDict_.lookup("contactGapTol"))),
    contactIterNum_(0),
    urDisp_
    (
        readScalar(normalContactModelDict_.lookup("relaxationFactor"))
    ),
    urTrac_
    (
        normalContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactorTraction",
            urDisp_
        )
    ),
    distanceMethod_(normalContactModelDict_.lookup("distanceMethod")),
    curTimeIndex_(-1),
    iCorr_(0),
    //slaveDispPrevIter_(slaveDisp_.size(), vector::zero),
    oscillationCorr_(normalContactModelDict_.lookup("oscillationCorrection")),
    smoothingSteps_(readInt(normalContactModelDict_.lookup("smoothingSteps"))),
    infoFreq_(readInt(normalContactModelDict_.lookup("infoFrequency"))),
    //debugFilePtr_(NULL),
    contactFilePtr_(NULL)
{
    if (dict.found("oldSlavePressure"))
    {
        Info<< "    Reading oldSlavePressure allowing restart of case" << endl;
        oldSlavePressure_ =
            vectorField
            (
                "oldSlavePressure",
                dict,
                mesh_.boundaryMesh()[slavePatchID].size()
            );

        // Add to dictionary for writing
        normalContactModelDict_.add("oldSlavePressure", oldSlavePressure_);
        // normalContactModelDict_.add("oldSlaveValueFrac", oldSlaveValueFrac_);
    }
    if (dict.found("oldSlaveValueFrac"))
    {
        Info<< "    Reading oldSlaveValueFrac allowing restart of case" << endl;
        oldSlaveValueFrac_ =
            symmTensorField
            (
                "oldSlaveValueFrac",
                dict,
                mesh_.boundaryMesh()[slavePatchID].size()
            );
    }

    // master proc open contact info file
    if (Pstream::master())
    {
        word masterName = mesh_.boundary()[masterPatchID].name();
        word slaveName = mesh_.boundary()[slavePatchID].name();
        contactFilePtr_ =
            new OFstream
            (
                fileName("normalContact_"+masterName+"_"+slaveName+".txt")
            );
        OFstream& contactFile = *contactFilePtr_;
        int width = 20;
        contactFile
            << "time";
        contactFile.width(width);
        contactFile
            << "iter";
        contactFile.width(width);
        contactFile
            << "contactFaces";
        contactFile.width(width);
        contactFile
            << "settledFaces";
        contactFile.width(width);
        contactFile
            << "minPen";
        contactFile.width(width);
        contactFile
            << "maxPress";
        contactFile.width(width);
        contactFile
            << "corrPoints" << endl;
    }

    // calc point points for missed vertices
    if (correctMissedVertices_)
    {
        calcSlavePointPoints(slavePatchID);
    }
}


dirichletNeumann::dirichletNeumann(const dirichletNeumann& nm)
:
    normalContactModel(nm),
    normalContactModelDict_(nm.normalContactModelDict_),
    mesh_(nm.mesh_),
    slaveDisp_(nm.slaveDisp_),
    slavePressure_(nm.slavePressure_),
    oldSlavePressure_(nm.oldSlavePressure_),
    zeroSlavePressure_(nm.zeroSlavePressure_),
    residual_(nm.residual_),
    maxNorm_(nm.maxNorm_),
    averagePenetration_(nm.averagePenetration_),
    nonLinearType_(nm.nonLinearType_),
    areaInContact_(nm.areaInContact_),
    slaveValueFrac_(nm.slaveValueFrac_),
    oldSlaveValueFrac_(nm.oldSlaveValueFrac_),
    limitPenetration_(nm.limitPenetration_),
    penetrationLimit_(nm.penetrationLimit_),
    limitPressure_(nm.limitPressure_),
    pressureLimit_(nm.pressureLimit_),
    correctMissedVertices_(nm.correctMissedVertices_),
    slavePointPointsPtr_(NULL),
    contactGapTol_(nm.contactGapTol_),
    contactIterNum_(nm.contactIterNum_),
    urDisp_(nm.urDisp_),
    urTrac_(nm.urTrac_),
    distanceMethod_(nm.distanceMethod_),
    curTimeIndex_(nm.curTimeIndex_),
    iCorr_(nm.iCorr_),
    oscillationCorr_(nm.oscillationCorr_),
    smoothingSteps_(nm.smoothingSteps_),
    infoFreq_(nm.infoFreq_),
    contactFilePtr_(NULL)
{
    if (nm.slavePointPointsPtr_)
    {
        slavePointPointsPtr_ = new labelListList(*nm.slavePointPointsPtr_);
    }

    if (nm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*nm.contactFilePtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

void dirichletNeumann::correct
(
    const primitiveFacePatch& masterFaceZonePatch,
    const primitiveFacePatch& slaveFaceZonePatch,
    const intersection::algorithm alg,
    const intersection::direction dir,
    word fieldName,
    const Switch orthotropic,
    const nonLinearGeometry::nonLinearType nonLinear,
    vectorField& slaveFaceNormals,
    const extendedGgiZoneInterpolation& zoneToZone
)
{
    const fvMesh& mesh = mesh_;
    iCorr_++;

    if
    (
        curTimeIndex_ != mesh.time().timeIndex()
        || nonLinearType_ !=
        word(mesh.solutionDict().subDict("solidMechanics").lookup("nonLinear"))
    )
    {
        curTimeIndex_ = mesh.time().timeIndex();
        iCorr_ = 0;
        maxNorm_ = SMALL;
        nonLinearType_ =
            word
            (
                mesh.solutionDict().subDict("solidMechanics").lookup("nonLinear")
            );
    }

    {
        //Info << "Correcting contact..." << flush;
        const label slavePatchIndex = slavePatchID();
        //const label masterPatchIndex = masterPatchID();
        contactIterNum_++;

        // calculate penetration distances from all slave faces
        // using either the point distance interpolated to the faces
        // or directly using the face distances

        //- Calculate intersection distances
        scalarField slaveDispMag(mesh.boundary()[slavePatchIndex].size(), 0.0);
        scalar minSlavePointPenetration = 0.0;
        scalarField globalSlavePointPenetration
            (
                slaveFaceZonePatch.points().size(), 0.0
            );
        scalarField slavePointPenetration
            (
                mesh.boundaryMesh()[slavePatchID()].nPoints(), 0.0
            );
        label numCorrectedPoints = 0;

        if (distanceMethod_ == "pointGGI")
        {
            // New distance method 'pointGGI" is fastest
            const scalarField& globalSlavePointPenetration =
                zoneToZone.slavePointDistanceToIntersection();

            // get local point values from global values
            forAll(slavePointPenetration, pointI)
            {
                // the local point values seem to be kept at the start of the
                // global field
                slavePointPenetration[pointI] =
                    globalSlavePointPenetration[pointI];

                //- when the master surface surrounds the slave (like the pelvis
                // and femur head) then the slave penetration can sometimes
                // calculate the distance through the femur head
                // to the pelvis which is wrong so I limit slavePenetration here
                if (limitPenetration_)
                {
                    if (slavePointPenetration[pointI] < penetrationLimit_)
                    {
                        slavePointPenetration[pointI] = 0.0;
                    }
                }
            }

            int num = 0;
            forAll(globalSlavePointPenetration, pointi)
            {
                if (globalSlavePointPenetration[pointi] < SMALL)
                {
                    num++;
                    averagePenetration_ += globalSlavePointPenetration[pointi];
                }
            }
            if (num > 0)
            {
                averagePenetration_ /= num;
            }

            // interpolate point distances to faces
            // set all positive penetrations to zero before interpolation
            primitivePatchInterpolation localSlaveInterpolator
                (
                    mesh.boundaryMesh()[slavePatchIndex]
                );
            slaveDispMag =
                localSlaveInterpolator.pointToFaceInterpolate<scalar>
                (
                    min(slavePointPenetration, 0.0)
                );

            // for visualisation
            slaveContactPointGap() = slavePointPenetration;
        }
        else if (distanceMethod_ == "point")
        {
            // create master to slave interpolation for calculation of distances
            PatchToPatchInterpolation<primitiveFacePatch,primitiveFacePatch>
                masterToSlavePatchToPatchInterpolator
                (
                    masterFaceZonePatch, // from zone
                    slaveFaceZonePatch, // to zone
                    alg,
                    dir
                );

            globalSlavePointPenetration =
                masterToSlavePatchToPatchInterpolator
                    .pointDistanceToIntersection();

            // get local point values from global values
            forAll(slavePointPenetration, pointI)
            {
                // the local point values seem to be kept at the start of the
                // global field
                slavePointPenetration[pointI] =
                    globalSlavePointPenetration[pointI];

                //- when the master surface surrounds the slave (like the pelvis
                // and femur head) then the slave penetration can sometimes
                // calculate the distance through the femur head
                // to the pelvis which is wrong so I limit slavePenetration here
                if (limitPenetration_)
                {
                    if (slavePointPenetration[pointI] < penetrationLimit_)
                    {
                        slavePointPenetration[pointI] = 0.0;
                    }
                }
            }

            if (correctMissedVertices_)
            {
                scalarField& slavePointPenetration_ = slavePointPenetration;
#               include "iterativePenaltyCorrectMissedVertices.H"
            }

            minSlavePointPenetration = gMin(globalSlavePointPenetration);

            // interpolate point distances to faces
            // set all positive penetrations to zero before interpolation
            primitivePatchInterpolation localSlaveInterpolator
                (
                    mesh.boundaryMesh()[slavePatchIndex]
                );
            slaveDispMag =
                localSlaveInterpolator.pointToFaceInterpolate<scalar>
                (
                    min(slavePointPenetration, 0.0)
                );

            // for visualisation
            slaveContactPointGap() = slavePointPenetration;
        }
        else if (distanceMethod_ == "face")
        {
            // create master to slave interpolation for calculation of distances
            PatchToPatchInterpolation<primitiveFacePatch, primitiveFacePatch>
                masterToSlavePatchToPatchInterpolator
                (
                    masterFaceZonePatch, // from zone
                    slaveFaceZonePatch, // to zone
                    alg,
                    dir
                );

            scalarField globalSlavePenetration =
                masterToSlavePatchToPatchInterpolator
                    .faceDistanceToIntersection();

            // get local point values from global values
            const label slavePatchStart
                = mesh.boundaryMesh()[slavePatchIndex].start();
            forAll(slaveDispMag, facei)
            {
                slaveDispMag[facei] =
                    globalSlavePenetration
                    [
                        mesh.faceZones()[slaveFaceZoneID()].whichFace
                        (
                            slavePatchStart + facei
                        )
                    ];

                // limit slavePenetration here
                if (limitPenetration_)
                {
                    if (slaveDispMag[facei] < penetrationLimit_)
                    {
                        slaveDispMag[facei] = 0.0;
                    }
                }
            }

            // we need the point distance to calculate the touch fraction
            // so we interpolate the face distances to points
            primitivePatchInterpolation localSlaveInterpolator
                (
                    mesh.boundaryMesh()[slavePatchIndex]
                );
            slavePointPenetration =
                localSlaveInterpolator.faceToPointInterpolate<scalar>
                (
                    slaveDispMag
                );

            // for visualisation
            slaveContactPointGap() = slavePointPenetration;

            // set all positive penetrations to zero
            slaveDispMag = min(slaveDispMag,0.0);

            minSlavePointPenetration = min(slavePointPenetration);
            reduce(minSlavePointPenetration, minOp<scalar>());
        }
        else
        {
            FatalError
                << "distanceMethod " << distanceMethod_ << " is unknown,\n"
                << "distanceMethod options are:\n    pointGGI\n    face"
                << nl << "    point"
                << abort(FatalError);
        }

        // calculate local deformed point and face normals
        // vectorField slaveFaceNormals(slaveDisp_.size(), vector::zero);
        //vectorField slavePointNormals
        //(slavePointPenetration.size(), vector::zero);

        // calculate slaveFaceNormals
        // we have a number of options:
        // we can use the slave deformed/undeformed normals
        // or the master deformed/undeformed normals interpolated to the slave
        // for large strain, we use the deformed normals
        // for small strain, we use the undeformed (or maybe deformed) normals
        //#       include "dirichletNeumannCalculateSlaveFaceNormals.H"
        calculateSlaveFaceNormals
        (
            slaveFaceNormals,
            nonLinear,
            masterFaceZonePatch,
            slaveFaceZonePatch,
            zoneToZone
        );

        // Calculate area in contact for local slave patch
        const faceList& slavePatchLocalFaces =
            mesh.boundaryMesh()[slavePatchIndex].localFaces();
        const pointField& slavePatchLocalPoints =
            mesh.boundaryMesh()[slavePatchIndex].localPoints();
        forAll (slavePatchLocalFaces, facei)
        {
            areaInContact_[facei] =
                slavePatchLocalFaces[facei].areaInContact
                (
                    slavePatchLocalPoints,
                    slavePointPenetration
                );
        }

        // set slave value fraction
        // fix normals of any face with non-zero increment of displacemet
        // and set displacement of faces inside contactGapTol to zero
        // count faces in contact
        int numSlaveContactFaces = 0;
        int numSlaveSettledFaces = 0;

        // const volVectorField& oldSlaveDispField =
        //     mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
        // const vectorField& oldSlaveDisp =
        //     oldSlaveDispField.boundaryField()[slavePatchIndex];
        //const scalarField oldSlaveDispMag = slaveFaceNormals & oldSlaveDisp;

        // set displacement increment
        forAll(areaInContact_, facei)
        {
            if (areaInContact_[facei] > SMALL)
            //if (areaInContact_[facei] > 1e-6)
            {
                numSlaveContactFaces++;

                // the edge of the contact area is often where convergence
                // difficulties occur.
                // using the touchFrac makes it difficult to converge.
                // the two versions give different convergence
                //slaveValueFrac_[facei] =
                //    areaInContact_[facei]*sqr(slaveFaceNormals[facei]);
                slaveValueFrac_[facei] = sqr(slaveFaceNormals[facei]);

                // Add contact tol
                slaveDispMag[facei] += contactGapTol_;

                if (slaveDispMag[facei] > -contactGapTol_)
                {
                    // count faces within gap tolerance
                    numSlaveSettledFaces++;
                    // apply extra relaxation
                    //slaveDispMag[facei] *= relaxFactor_;
                    //slaveDispMag[facei] = 0.0;
                }
            }
            else
            {
                // face not in contact
                slaveValueFrac_[facei] = symmTensor::zero;
                //slaveDispMag[facei] = max(0.0, slaveDispMag[facei]);
                //slaveDispMag[facei] = 0.0;

                // no relaxation on pressure or slave disp
                //oldSlavePressure_[facei] = vector::zero;
                //slaveDisp_[facei] = vector::zero;
            }
        }

        reduce(numSlaveContactFaces, sumOp<int>());
        reduce(numSlaveSettledFaces, sumOp<int>());

        // fixed under-relaxation
        //slaveDispMag *= relaxFactor_; // WAS THIS

        // set slave displacement increment
        //slaveDisp_ = slaveFaceNormals*slaveDispMag; // WAS THIS
        // slaveDisp_ =
        //     relaxFactor_*slaveFaceNormals*slaveDispMag
        //     + (1.0 - relaxFactor_)*slaveDisp_;
        slaveDisp_ += urDisp_*slaveFaceNormals*slaveDispMag;

        // slaveDisp is a correction to the current
        // patch U or DU so we add it to the current U/DU
        // const volVectorField& oldSlaveDispField =
        //   mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
        // const vectorField& oldSlaveDisp =
        //oldSlaveDispField.boundaryField()[slavePatchIndex];
        //slaveDisp_ += oldSlaveDisp; // WAS THIS

        // Remove tengential component
        slaveDisp_ = slaveFaceNormals*(slaveFaceNormals & slaveDisp_);

        //--------Try to remove any oscillations----------//
        if (oscillationCorr_)
        {
            if (smoothingSteps_ < 1)
            {
                FatalError
                    << "SmoothingSteps must be greater than or equal to 1"
                    << abort(FatalError);
            }

            // interpolate face values to points then interpolate back
            // this essentially smooths the field
            primitivePatchInterpolation localSlaveInterpolator
                (
                    mesh.boundaryMesh()[slavePatchIndex]
                );
            vectorField slaveDispPoints
                (
                    mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero
                );

            for (int i=0; i<smoothingSteps_; i++)
            {
                slaveDispPoints =
                    localSlaveInterpolator.faceToPointInterpolate<vector>
                    (
                        slaveDisp_
                    );
                slaveDisp_ =
                    localSlaveInterpolator.pointToFaceInterpolate<vector>
                    (
                        slaveDispPoints
                    );

                // make sure no tangential component
                slaveDisp_ = slaveFaceNormals*(slaveFaceNormals & slaveDisp_);
            }
        }

        // Calculate current slave traction
        {
            const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
            const fvPatchField<tensor>& gradField =
                slavePatch.lookupPatchField<volTensorField, tensor>
                (
                    "grad(" + fieldName + ")"
                );

            slavePressure_ =
                tractionBoundaryGradient().traction
                (
                    gradField,                 // grad field
                    fieldName,                 // working field name
                    "U",                       // total field name
                    slavePatch,                // polyPatch
                    bool(fieldName == "DU")    // incremental
                );

            // Set traction to zero on faces not in contact
            // and set tensile stresses to zero
            //const scalar maxMagSlavePressure = gMax(mag(slavePressure_));
            scalar maxMagSlavePressure = 0.0;
            if (slavePressure_.size() > 0)
            {
                maxMagSlavePressure = max(mag(slavePressure_));
            }
            reduce(maxMagSlavePressure, maxOp<scalar>());

            forAll(areaInContact_, facei)
            {
                if
                (
                    // areaInContact_[facei] < SMALL
                    // || ( (slaveFaceNormals[facei] & slavePressure_[facei])
                    //     > 1e-3*maxMagSlavePressure)
                    (slaveFaceNormals[facei] & slavePressure_[facei]) > SMALL
                )
                {
                    slavePressure_[facei] = vector::zero;
                    //oldSlavePressure_[facei] = vector::zero;
                    slaveValueFrac_[facei] = symmTensor::zero;
                    //oldSlaveValueFrac_[facei] = symmTensor::zero;
                    // slaveDisp_[facei] =
                    //     slaveFaceNormals[facei]
                    //     *(slaveFaceNormals[facei] & oldSlaveDisp[facei]);
                }
            }

            // Under-relax traction
            // slavePressure_ =
            //      relaxFactor_*slavePressure_
            //      + (1-relaxFactor_)*oldSlavePressure_;
            slavePressure_ =
                urTrac_*slavePressure_ + (1 - urTrac_)*oldSlavePressure_;
            oldSlavePressure_ = slavePressure_;

            // Remove any shear tractions
            slavePressure_ =
                slaveFaceNormals*(slaveFaceNormals & slavePressure_);

            // Limiting pressure to can help convergence if max pressure is
            // known
            if (limitPressure_)
            {
                forAll(slavePressure_, facei)
                {
                    if
                    (
                        (slaveFaceNormals[facei] & slavePressure_[facei])
                        < -pressureLimit_
                    )
                    {
                        slavePressure_[facei] =
                            -pressureLimit_*slaveFaceNormals[facei];
                    }
                }
            }

            // Calculate residual
            //maxNorm_ = max(gSum(magSqr(slavePressure_)), maxNorm_);
            //residual_ =
            //    gSum(magSqr(slavePressure_ - oldSlavePressure_))/maxNorm_;

            residual_ =
                 gSum(magSqr(slavePressure_ - oldSlavePressure_))
                 /(gSum(magSqr(oldSlavePressure_)) + SMALL);
            //maxNorm_ = max(residual_, maxNorm_);
        }


        // In parallel, the log is poluted with warnings that
        // I am getting max of a list of size zero so
        // I will get the max of procs which have some
        // of the slave faces
        //scalar maxMagMasterTraction = gMax(mag(slavePressure_))
        scalar maxMagMasterTraction = 0.0;
        if (slavePressure_.size() > 0)
        {
            maxMagMasterTraction = max(mag(slavePressure_));
        }
        reduce(maxMagMasterTraction, maxOp<scalar>());

        // Under-relax value fraction - disabled
        slaveValueFrac_ =
              urDisp_*slaveValueFrac_
              + (1.0 - urDisp_)*oldSlaveValueFrac_;
        oldSlaveValueFrac_ = slaveValueFrac_;

        if (Pstream::master() && (contactIterNum_ %  infoFreq_ == 0))
        {
            OFstream& contactFile = *contactFilePtr_;
            int width = 20;
            contactFile
                << mesh.time().value();
            contactFile.width(width);
            contactFile
                << contactIterNum_;
            contactFile.width(width);
            contactFile
                << numSlaveContactFaces;
            contactFile.width(width);
            contactFile
                << numSlaveSettledFaces;
            contactFile.width(width);
            contactFile
                << minSlavePointPenetration;
            contactFile.width(width);
            contactFile
                << maxMagMasterTraction;
            contactFile.width(width);
            contactFile
                << numCorrectedPoints;
            contactFile << endl;
        }
    }
}

void dirichletNeumann::calcSlavePointPoints(const label slavePatchID)
{
    // calculate patch pointPoints
    // code adapted from enrichedPatch calcPointPoints

    // Calculate point-point addressing
    if (slavePointPointsPtr_)
    {
        FatalErrorIn("void contactPatchPair::calcSlavePointPoints() const")
            << "Point-point addressing already calculated."
                << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    // Algorithm:
    // Go through all faces and add the previous and next point as the
    // neighbour for each point. While inserting points, reject the
    // duplicates (as every internal edge will be visited twice).
    List<DynamicList<label, primitiveMesh::edgesPerPoint_> >
        pp(mesh.boundaryMesh()[slavePatchID].meshPoints().size());

    const faceList& lf = mesh.boundaryMesh()[slavePatchID].localFaces();

    register bool found = false;

    forAll (lf, faceI)
    {
        const face& curFace = lf[faceI];

        forAll (curFace, pointI)
        {
            DynamicList<label, primitiveMesh::edgesPerPoint_>&
                curPp = pp[curFace[pointI]];

            // Do next label
            label next = curFace.nextLabel(pointI);

            found = false;

            forAll (curPp, i)
            {
                if (curPp[i] == next)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                curPp.append(next);
            }

            // Do previous label
            label prev = curFace.prevLabel(pointI);
            found = false;

            forAll (curPp, i)
            {
                if (curPp[i] == prev)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                curPp.append(prev);
            }
        }
    }

    // Re-pack the list
    slavePointPointsPtr_ = new labelListList(pp.size());
    labelListList& ppAddr = *slavePointPointsPtr_;

    forAll (pp, pointI)
    {
        ppAddr[pointI].transfer(pp[pointI].shrink());
    }
}


scalar dirichletNeumann::residual()
{
    // if (averagePenetration_ > 5*contactGapTol_)
    // {
    //     return 1.0;
    // }

    return residual_;
};

void dirichletNeumann::writeDict(Ostream& os) const
{
    word keyword(name()+"NormalModelDict");

    os.writeKeyword(keyword)
        << normalContactModelDict_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
