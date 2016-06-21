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
    iterativePenalty

\*---------------------------------------------------------------------------*/

#include "iterativePenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(iterativePenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, iterativePenalty, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

iterativePenalty::iterativePenalty
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID,
    const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
    const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
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
  slaveValueFrac_(mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero),
  globalSlavePointPenetration_(slaveFaceZonePatch.points().size(), 0.0),
  slavePointPenetration_
  (mesh_.boundaryMesh()[slavePatchID].localPoints().size(), 0.0),
  limitPenetration_(normalContactModelDict_.lookup("limitPenetration")),
  penetrationLimit_
  (readScalar(normalContactModelDict_.lookup("penetrationLimit"))),
  correctMissedVertices_
  (normalContactModelDict_.lookup("correctMissedVertices")),
  slavePointPointsPtr_(NULL),
  penaltyFactorPtr_(NULL),
  returnPenaltyFactorPtr_(NULL),
  penaltyScale_(readScalar(normalContactModelDict_.lookup("penaltyScale"))),
  returnScale_(readScalar(normalContactModelDict_.lookup("returnScale"))),
  totalSlavePointTrac_
  (
   IOobject
   (
    "totalSlavePointTraction",
    mesh_.time().timeName(),
    mesh_,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   vectorField(slavePointPenetration_.size(), vector::zero)
   ),
  contactGapTol_(readScalar(normalContactModelDict_.lookup("contactGapTol"))),
  contactIterNum_(0),
  oscillationCorr_(normalContactModelDict_.lookup("oscillationCorrection")),
  oscillationCorrFac_
  (readScalar(normalContactModelDict_.lookup("oscillationCorrectionFactor"))),
  contactFilePtr_(NULL)
{
  // master proc open contact info file
  if (Pstream::master())
    {
      word masterName = mesh_.boundary()[masterPatchID].name();
      word slaveName = mesh_.boundary()[slavePatchID].name();
      contactFilePtr_ =
          new OFstream
          (fileName("normalContact_"+masterName+"_"+slaveName+".txt"));
      OFstream& contactFile = *contactFilePtr_;
      int width = 20;
      contactFile << "Time";
      contactFile.width(width);
      contactFile << "ContactIter";
      contactFile.width(width);
      contactFile << "PenaltyScale";
      contactFile.width(width);
      contactFile << "slaveContactVer";
      contactFile.width(width);
      contactFile << "penetration";
      contactFile.width(width);
      contactFile << "slaveVerReduce";
      contactFile.width(width);
      contactFile << "maxSlaveTrac" << endl;
    }

  // calc point points for missed vertices
  if (correctMissedVertices_)
    {
      calcSlavePointPoints(slavePatchID);
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

  void iterativePenalty::correct
  (
   // const PrimitivePatchInterpolation
   // < PrimitivePatch<face, List, pointField> >& slaveInterpolator,
   const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
   const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
   const intersection::algorithm alg,
   const intersection::direction dir,
   word fieldName,
   //const label slaveFaceZoneID,
   const Switch orthotropic,
   const word nonLinear,
   vectorField& slaveFaceNormals,
   GGIInterpolation< PrimitivePatch< face, List, pointField >,
             PrimitivePatch< face, List, pointField >
             >* ggiInterpolatorPtr
   )
  {
    //    Info << "Correcting contact..." << flush;

    //---------------------PRELIMINARIES---------------------------------//
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    //const label masterPatchIndex = masterPatchID();
    scalar maxMagSlaveTraction = 0.0;
    contactIterNum_++;


    //--------------------CALCULATE SLAVE PENETRATION----------------------//
    // create master to slave interpolation for calculation of distances
   PatchToPatchInterpolation<
     PrimitivePatch<
         face, List, pointField
         >, PrimitivePatch<face, List, pointField>
       > masterToSlavePatchToPatchInterpolator
       (
    masterFaceZonePatch, // from zone
    slaveFaceZonePatch, // to zone
    alg,
    dir
    );
  //- calculate intersection distances
  //- this is the slowest part of the contact correction especially
  // when the slavePatch
  //- has many points. parallelisation of this step should be considered.
   globalSlavePointPenetration_
    = masterToSlavePatchToPatchInterpolator.pointDistanceToIntersection();

  forAll(slavePointPenetration_, pointI)
    {
      // the local point values seem to be kept at the start of the global field
      slavePointPenetration_[pointI] =
    globalSlavePointPenetration_
    [
     pointI
     ];

      //- when the master surface surrounds the slave (like the pelvis and
      // femur head) then
      //- the slave penetration can sometimes calculate the distance through
      // the femur head
      //- to the pelvis which is wrong so I limit slavePenetration here
      //- i should add a limitPenetration switch here if (limitPenetration)
      if (limitPenetration_)
    {
      if (slavePointPenetration_[pointI] < penetrationLimit_)
        {
          slavePointPenetration_[pointI] = 0.0;
        }
    }
    }

  label numCorrectedPoints = 0;
  if (correctMissedVertices_)
    {
#     include "iterativePenaltyCorrectMissedVertices.H"
    }

   // calculate slaveFaceNormals
   // we have a number of options:
   // we can use the slave deformed/undeformed normals
   // or the master deformed/undeformed normals interpolated to the slave
   // for large strain, we use the deformed normals
   // for small strain, we use the undeformed (or maybe deformed) normals
#  include "dirichletNeumannCalculateSlaveFaceNormals.H"


  //------CALCULATE SLAVE VERTEX TRACTIONS BASED ON PENETRATION-------------//
  // use the deformed normals from the faceZonePatches
  const vectorField& globalSlavePointNormals =
      slaveFaceZonePatch.pointNormals();
  vectorField slavePointNormals(slavePointPenetration_.size(), vector::zero);
  forAll(slavePointNormals, pointi)
    {
      slavePointNormals[pointi] = globalSlavePointNormals[pointi];
    }

  int numSlaveContactPoints = 0;
  int numSlaveContactPointsReducing = 0;
  int numSlavesUpdated = 0;

  //- so the procs know the global min
  scalar minSlavePointPenetration = gMin(globalSlavePointPenetration_);
  scalar penaltyFac = penaltyFactor();
  scalar returnPenaltyFac = returnPenaltyFactor();

  // calculate traction increments
  forAll(totalSlavePointTrac_, pointI)
    {
      // if a point has penetrated (i.e. if the penetration is negative),
      // add an increment of traction to it relative to the penetration
      if (slavePointPenetration_[pointI] < -contactGapTol_)
    {
      //contactStep = true;
      numSlaveContactPoints++; // count points in contact
      numSlavesUpdated++;
      //- traction is linearly dependent on penetration
      totalSlavePointTrac_[pointI] +=
        ( slavePointNormals[pointI] * penaltyFac
          * slavePointPenetration_[pointI] );
    }
      else if (slavePointPenetration_[pointI] < 0.0)
    {
      numSlaveContactPoints++; // count points in contact
    }
      // else if penetration is positive and there is a positive
      // pressure (negative traction) still
      // on the point, then slowly reduce the pressure
      else if ((totalSlavePointTrac_[pointI] & slavePointNormals[pointI]) < 0.0)
    {
      numSlavesUpdated++;
      numSlaveContactPointsReducing++;
      // point tractions must be reduced slowly

      totalSlavePointTrac_[pointI] +=
        ( slavePointNormals[pointI] * returnPenaltyFac
          *slavePointPenetration_[pointI] );

      // if a tensile traction develops
      if ((totalSlavePointTrac_[pointI] & slavePointNormals[pointI]) > 0.0)
        {
          totalSlavePointTrac_[pointI] = vector::zero;
        }
    }
    }


  //--------INTERPOLATE SLAVE POINT TRACTION TO SLAVE FACES----------//
  //- create local patch interpolation
  //- no need to interpolate using the entire slave face zone patch
  primitivePatchInterpolation localSlaveInterpolator
      (mesh.boundaryMesh()[slavePatchIndex]);
  slavePressure_ =
    localSlaveInterpolator.pointToFaceInterpolate<vector>
    (
     totalSlavePointTrac_
     );

  //--------Try to remove any oscillations----------//
   if (oscillationCorr_)
   {
     correctOscillations(slavePressure_, slaveFaceZonePatch);
   }

   // in parallel, the log is poluted with warnings that
   // I am getting max of a list of size zero so
   // I will get the max of procs which have some
   // of the slave faces
   // maxMagSlaveTraction = gMax(mag(slavePressure_));
  if (slavePressure_.size() > 0) maxMagSlaveTraction = max(mag(slavePressure_));
  reduce(maxMagSlaveTraction, maxOp<scalar>());


  //--------MASTER PROCS WRITES CONTACT INFO FILE----------//
  reduce(numSlaveContactPoints, sumOp<int>());
  reduce(numSlaveContactPointsReducing, sumOp<int>());

  if (Pstream::master())
    {
      OFstream& contactFile = *contactFilePtr_;
      int width = 15;
      contactFile << mesh.time().value();
      contactFile.width(width);
      contactFile << contactIterNum_;
      contactFile.width(width);
      contactFile << penaltyScale_;
      contactFile.width(width);
      contactFile << numSlaveContactPoints;
      contactFile.width(width);
      contactFile << minSlavePointPenetration;
      contactFile.width(width);
      contactFile << numSlaveContactPointsReducing;
      contactFile.width(width);
      contactFile << maxMagSlaveTraction;
      contactFile.width(width);
      contactFile << " " << numCorrectedPoints << endl;
    }
  //  Info << "\tdone" << endl;
  }


  void iterativePenalty::correctOscillations
  (
   vectorField& slavePressure,
   const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
   )
  {
    // oscillations sometimes appear in normal contact pressure
    // so we will try to limit them here
    // we will weight the current slavePressure the average of the
    // neighbours using the weight oscillationCorrectionFactor_

    //Pout << "Applying contact oscillation correction..." << flush;

    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    const labelListList& faceFaces = slaveFaceZonePatch.faceFaces();

    // create global slavePressure
    vectorField globalSlavePressure(slaveFaceZonePatch.size(), vector::zero);
    const label slavePatchStart
      = mesh.boundaryMesh()[slavePatchIndex].start();
    forAll(slavePressure, i)
      {
          globalSlavePressure
              [
                  mesh.faceZones()[slaveFaceZoneID()
                      ].whichFace(slavePatchStart + i)] =
      slavePressure[i];
      }
    // sum because each face is only on one proc
    reduce(globalSlavePressure, sumOp<vectorField>());

    // smooth slavePressure with face face disps
    forAll(faceFaces, facei)
      {
    //label numFaceFaces = faceFaces[facei].size();

    if (mag(globalSlavePressure[facei]) > SMALL)
      //&&
      //numFaceFaces > 1) // don't smooth end/corner faces
      {
        vector avPress = vector::zero;
        int numNei = 0;
        forAll(faceFaces[facei], ffi)
          {
        label faceFace = faceFaces[facei][ffi];

        avPress += globalSlavePressure[faceFace];
        numNei++;
        }

        // if (numNei == 0)
        //   {
        //      avPress = globalSlavePressure[facei];
        //   }
        // else
        //   {
        //      avPress /= faceFaces[facei].size();
        //   }
        avPress /= numNei;

        // if (numFaceFaces == 1)
        //   {
        //      // for corner/end faces, decrease the weight of the neighbours
        //      avPress += globalSlavePressure[facei];
        //      avPress /= 2;
        //   }

        // weighted-average with face-faces
        globalSlavePressure[facei] =
          oscillationCorrFac_*globalSlavePressure[facei]
            + (1.0-oscillationCorrFac_)*avPress;
      }
      }

    // convert global back to local
    forAll(slavePressure, facei)
       {
     slavePressure[facei] =
       globalSlavePressure
       [
        mesh.faceZones()[slaveFaceZoneID()].whichFace(slavePatchStart + facei)
        ];
       }

    //Pout << "\tdone" << endl;
  }


  void iterativePenalty::calcSlavePointPoints(const label slavePatchID)
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


  void iterativePenalty::calcPenaltyFactor()
  {
    // set penalty factor and return factor
    // approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale
    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();
    const constitutiveModel& rheology =
      mesh_.objectRegistry::lookupObject<constitutiveModel>
        ("rheologyProperties");
    scalarField masterMu =
      rheology.mu()().boundaryField()[masterPatchIndex];
    scalarField slaveMu =
      rheology.mu()().boundaryField()[slavePatchIndex];
    scalarField masterLambda =
      rheology.lambda()().boundaryField()[masterPatchIndex];
    scalarField slaveLambda =
      rheology.lambda()().boundaryField()[slavePatchIndex];

    // avarage contact patch bulk modulus
    scalar masterK = gAverage(masterLambda + (2.0/3.0)*masterMu);
    scalar slaveK = gAverage(slaveLambda + (2.0/3.0)*slaveMu);
    scalar bulkModulus = 0.5*(masterK+slaveK);

    // average contact patch face area
    scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    scalar faceArea = 0.5*(masterMagSf+slaveMagSf);

    // average contact patch cell volume
    scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    const volScalarField::DimensionedInternalField & V = mesh_.V();
    {
      const unallocLabelList& faceCells =
          mesh_.boundary()[masterPatchIndex].faceCells();
      forAll(mesh_.boundary()[masterPatchIndex], facei)
    {
      masterV[facei] = V[faceCells[facei]];
    }
    }
    {
      const unallocLabelList& faceCells =
          mesh_.boundary()[slavePatchIndex].faceCells();
      forAll(mesh_.boundary()[slavePatchIndex], facei)
    {
      slaveV[facei] = V[faceCells[facei]];
    }
    }
    scalar cellVolume = 0.5*(gAverage(masterV)+gAverage(slaveV));

    // approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    penaltyFactorPtr_ =
        new scalar(penaltyScale_*bulkModulus*faceArea/cellVolume);
  }

  void iterativePenalty::writeDict(Ostream& os) const
  {
    word keyword(name()+"NormalModelDict");
    os.writeKeyword(keyword)
      << normalContactModelDict_;
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
