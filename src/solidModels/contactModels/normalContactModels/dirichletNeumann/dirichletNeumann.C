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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(dirichletNeumann, 0);
  addToRunTimeSelectionTable(normalContactModel, dirichletNeumann, dictionary);

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
    oldSlaveDispMag_(mesh_.boundaryMesh()[slavePatchID].size(), 0.0),
    slavePressure_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    oldSlavePressure_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    touchFraction_(mesh_.boundaryMesh()[slavePatchID].size(), 0.0),
    slaveValueFrac_(mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero),
    oldSlaveValueFrac_
    (
        mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero
    ),
    limitPenetration_(normalContactModelDict_.lookup("limitPenetration")),
    penetrationLimit_
    (
        readScalar(normalContactModelDict_.lookup("penetrationLimit"))
    ),
    limitPressure_(normalContactModelDict_.lookup("limitPressure")),
    pressureLimit_(readScalar(normalContactModelDict_.lookup("pressureLimit"))),
    settleContact_
    (
        normalContactModelDict_.lookupOrDefault<Switch>
        (
             "settleContact",
             false
        )
    ),
    settleIterationNumber_
    (
        normalContactModelDict_.lookupOrDefault<label>
        (
            "settleIterationNumber",
            1000
        )
    ),
    correctMissedVertices_
    (
        normalContactModelDict_.lookup("correctMissedVertices")
    ),
    slavePointPointsPtr_(NULL),
    contactGapTol_(readScalar(normalContactModelDict_.lookup("contactGapTol"))),
    contactIterNum_(0),
    relaxFactor_(readScalar(normalContactModelDict_.lookup("relaxationFactor"))),
    distanceMethod_(normalContactModelDict_.lookup("distanceMethod")),
    aitkenRelaxation_(normalContactModelDict_.lookup("aitkenRelaxation")),
    curTimeIndex_(-1),
    iCorr_(0),
    aitkenRes0_(1.0),
    aitkenTheta_(relaxFactor_),
    aitkenDelta_(slaveDisp_.size(), vector::zero),
    aitkenDeltaPrevIter_(slaveDisp_.size(), vector::zero),
    slaveDispPrevIter_(slaveDisp_.size(), vector::zero),
    oscillationCorr_(normalContactModelDict_.lookup("oscillationCorrection")),
    smoothingSteps_(readInt(normalContactModelDict_.lookup("smoothingSteps"))),
    infoFreq_(readInt(normalContactModelDict_.lookup("infoFrequency"))),
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
        contactFile << "time";
        contactFile.width(width);
        contactFile << "iter";
        contactFile.width(width);
        contactFile << "contactFaces";
        contactFile.width(width);
        contactFile << "settledFaces";
        contactFile.width(width);
        contactFile << "minPen";
        contactFile.width(width);
        contactFile << "maxPress";
        contactFile.width(width);
        contactFile << "corrPoints" << endl;
    }

    // calc point points for missed vertices
    if (correctMissedVertices_)
    {
        calcSlavePointPoints(slavePatchID);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void dirichletNeumann::correct
(
    const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
    const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
    const intersection::algorithm alg,
    const intersection::direction dir,
    word fieldName,
    const Switch orthotropic,
    const word nonLinear,
    vectorField& slaveFaceNormals,
    GGIInterpolation
    <
        PrimitivePatch< face, List, pointField >,
        PrimitivePatch< face, List, pointField >
    >* ggiInterpolatorPtr
)
{
    if (!settleContact_ || (iCorr_ < settleIterationNumber_))
    {

    //Info << "Correcting contact..." << flush;
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    //const label masterPatchIndex = masterPatchID();
    contactIterNum_++;

    // calculate penetration distances from all slave faces
    // using either the point distance interpolated to the faces
    // or directly using the face distances

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
   //- this is the slowest part of the contact correction
   // especially when the slavePatch
   //- has many points.
   // parallelisation of this step should be considered.
   scalarField slaveDispMag(mesh.boundary()[slavePatchIndex].size(), 0.0);
   scalar minSlavePointPenetration = 0.0;
   scalarField globalSlavePointPenetration
       (slaveFaceZonePatch.points().size(), 0.0);
   scalarField slavePointPenetration
       (mesh.boundaryMesh()[slavePatchID()].nPoints(), 0.0);
   label numCorrectedPoints = 0;
   if (distanceMethod_ == "point")
   {
       // pointDistanceToIntersection() sometimes gives a seg fault when the
       // momentum equation diverges

       globalSlavePointPenetration =
           masterToSlavePatchToPatchInterpolator.pointDistanceToIntersection();

       // get local point values from global values
       forAll(slavePointPenetration, pointI)
       {
           // the local point values seem to be kept at the start
           // of the global field
           slavePointPenetration[pointI] =
               globalSlavePointPenetration[pointI];

           //- when the master surface surrounds the slave (like the pelvis and
           // femur head) then
           //- the slave penetration can sometimes calculate the distance through
           // the femur head
           //- to the pelvis which is wrong so I limit slavePenetration here
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
#          include "iterativePenaltyCorrectMissedVertices.H"
       }

       minSlavePointPenetration = gMin(globalSlavePointPenetration);

       // interpolate point distances to faces
       // set all positive penetrations to zero before interpolation
       primitivePatchInterpolation localSlaveInterpolator
           (mesh.boundaryMesh()[slavePatchIndex]);
       slaveDispMag =
     localSlaveInterpolator.pointToFaceInterpolate<scalar>
     (
         min(slavePointPenetration, scalar(0))
     );

       // for visualisation
       slaveContactPointGap() = slavePointPenetration;
     }
   else if (distanceMethod_ == "face")
     {
       scalarField globalSlavePenetration
     = masterToSlavePatchToPatchInterpolator.faceDistanceToIntersection();

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

       //- when the master surface surrounds the slave (like the pelvis
       // and femur head) then
       //- the slave penetration can sometimes calculate the distance
       // through the femur head
       //- to the pelvis which is wrong so I limit slavePenetration here
       if (limitPenetration_)
         {
           if (slaveDispMag[facei] < penetrationLimit_)
         {
           slaveDispMag[facei] = 0.0;
         }
         }
     }

       // under-relax slaveDispMag
       slaveDispMag =
           relaxFactor_*slaveDispMag
           + (1.0 - relaxFactor_)*oldSlaveDispMag_;
       oldSlaveDispMag_ = slaveDispMag;

       // we need the point distance to calculate the touch fraction
       // so we interpolate the face distances to points
       primitivePatchInterpolation localSlaveInterpolator
           (mesh.boundaryMesh()[slavePatchIndex]);
       slavePointPenetration =
     localSlaveInterpolator.faceToPointInterpolate<scalar>
     (
      slaveDispMag
      );

       // for visualisation
       slaveContactPointGap() = slavePointPenetration;

       // set all positive penetrations to zero
       slaveDispMag = min(slaveDispMag, scalar(0));

       minSlavePointPenetration = min(slavePointPenetration);
       reduce(minSlavePointPenetration, minOp<scalar>());
     }
   else
     {
       FatalError << "distanceMethod " << distanceMethod_ << " is unknown,\n"
     "distanceMethod options are:\n\tface\n\tpoint"
          << exit(FatalError);
     }

   // calculate local deformed point and face normals
   // vectorField slaveFaceNormals(slaveDisp_.size(), vector::zero);
   //vectorField slavePointNormals(slavePointPenetration.size(), vector::zero);

   // calculate slaveFaceNormals
   // we have a number of options:
   // we can use the slave deformed/undeformed normals
   // or the master deformed/undeformed normals interpolated to the slave
   // for large strain, we use the deformed normals
   // for small strain, we use the undeformed (or maybe deformed) normals
#  include "dirichletNeumannCalculateSlaveFaceNormals.H"

   // calculate area in contact for local slave patch
   const faceList& slavePatchLocalFaces =
     mesh.boundaryMesh()[slavePatchIndex].localFaces();
   const pointField& slavePatchLocalPoints =
     mesh.boundaryMesh()[slavePatchIndex].localPoints();
   forAll (slavePatchLocalFaces, facei)
     {
       touchFraction_[facei] =
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

   const volVectorField& oldSlaveDispField =
     mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
   const vectorField& oldSlaveDisp =
       oldSlaveDispField.boundaryField()[slavePatchIndex];
   const scalarField oldSlaveDispMag = slaveFaceNormals & oldSlaveDisp;

   // set displacement increment
   forAll(touchFraction_, facei)
     {
       if (touchFraction_[facei] > SMALL)
     {
       numSlaveContactFaces++;

       // the edge of the contact area is often where convergence
       // difficulties occur.
       // using the touchFrac makes it difficult to converge.
       // the two versions give different convergence
       //slaveValueFrac_[facei] =
       //touchFraction_[facei]*sqr(slaveFaceNormals[facei]);
       slaveValueFrac_[facei] = sqr(slaveFaceNormals[facei]);

       slaveDispMag[facei] += contactGapTol_;

       // count faces within gap tolerance
       if (slaveDispMag[facei] > -contactGapTol_)
         //if (slaveDispMag[facei] > 0.0)
         {
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
       slaveDispMag[facei] = 0.0;
     }
     }

   reduce(numSlaveContactFaces, sumOp<int>());
   reduce(numSlaveSettledFaces, sumOp<int>());

   // if it is a new time step then reset iCorr
   iCorr_++;
   if (curTimeIndex_ != mesh.time().timeIndex())
     {
       curTimeIndex_ = mesh.time().timeIndex();
       iCorr_ = 0;
     }


   // under-relax the displacement increment
   if (aitkenRelaxation_)
     {
       // Aitken's method for under-relaxation
   if (curTimeIndex_ != mesh.time().timeIndex())
     {
       //Info << "aitkenRes0_ is " << aitkenRes0_ << endl;
       aitkenRes0_ = 1.0;
       aitkenTheta_ = relaxFactor_;
     }

       if (iCorr_ == 1)
     {
       //const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
       //const fvPatchField<vector>& dispField =
       //slavePatch.lookupPatchField<volVectorField, vector>(fieldName);
       const volVectorField& dispField =
         mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
       //all residuals will be normalised to the max mag of the DU/U field
       aitkenRes0_ = gMax(mag(dispField.internalField()));
       //Info << "updating aitkenRes0_ to " << aitkenRes0_ << endl;
     }

       // update delta
       // normalised with repsect of aitkenRes0 to avoid very small numbers
       slaveDisp_ = slaveFaceNormals*slaveDispMag;
       aitkenDelta_ = (slaveDisp_ - slaveDispPrevIter_) / aitkenRes0_;

       // first iteration set to relaxFactor
       if (iCorr_ > 0)
     {
       //Info << "slaveDisp_ " << slaveDisp_ << endl;
       vectorField b = aitkenDelta_ - aitkenDeltaPrevIter_;
       // scalar sumMagB = gSum(mag(b));
       scalar sumMagB = gSum(magSqr(b));
       if (sumMagB < SMALL)
         {
           Warning
               << "Dirichlet Normal Contatc Aitken under-relaxation: "
               << "denominator less then SMALL"
               << endl;
           sumMagB += SMALL;
         }
       aitkenTheta_ = -aitkenTheta_*
         gSum(aitkenDeltaPrevIter_ & b)
         /
         sumMagB;
     }

       // update disp with Aitken correction
       slaveDisp_ += aitkenTheta_*aitkenDelta_*aitkenRes0_;

       // update previous iteration value
       slaveDispPrevIter_ = slaveDisp_;
       aitkenDeltaPrevIter_ = aitkenDelta_;
     }
   else
     {
       // fixed under-relaxation
       slaveDispMag *= relaxFactor_;

       // set slave displacement increment
       slaveDisp_ = slaveFaceNormals*slaveDispMag;
     }



   // slaveDisp is a correction to the current
   // patch U or DU so we add it to the current U/DU
   // const volVectorField& oldSlaveDispField =
   //   mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
   // const vectorField& oldSlaveDisp =
   // oldSlaveDispField.boundaryField()[slavePatchIndex];
   slaveDisp_ += oldSlaveDisp;

   // remove tengential component
   slaveDisp_ = slaveFaceNormals*(slaveFaceNormals & slaveDisp_);

   //--------Try to remove any oscillations----------//
   if (oscillationCorr_)
   {
     if (smoothingSteps_ < 1)
       {
     FatalError << "smoothingSteps must be greater than or equal to 1"
            << exit(FatalError);
       }

     // interpolate face values to points then interpolate back
     // this essentially smooths the field
     primitivePatchInterpolation localSlaveInterpolator
         (mesh.boundaryMesh()[slavePatchIndex]);
     vectorField slaveDispPoints
         (mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero);

     for (int i=0; i<smoothingSteps_; i++)
     {
         slaveDispPoints =
             localSlaveInterpolator.faceToPointInterpolate<vector>(slaveDisp_);

         slaveDisp_ = localSlaveInterpolator.pointToFaceInterpolate<vector>
             (slaveDispPoints);

         // make sure no tangential component
         slaveDisp_ = slaveFaceNormals*(slaveFaceNormals & slaveDisp_);
     }
   }


   // write slaveDisps to file for debugging
   //    {
   //      OFstream& debugFile = *debugFilePtr_;
   //      forAll(slaveDisp_, facei)
   //        {
   //    if (facei > 35)
   //      {
   //        debugFile << mag(slaveDisp_[facei]) << " ";
   //      }
   //        }
   //      debugFile << endl;
   //    }


   // calculate current slave traction
   {
       const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
       const fvPatchField<tensor>& gradField =
           slavePatch.lookupPatchField<volTensorField, tensor>
           (
               "grad(" + fieldName + ")"
           );

       bool incremental(fieldName == "DU");

       slavePressure_ = tractionBoundaryGradient::traction
       (
           gradField,
           fieldName,
           "U",
           slavePatch,
           orthotropic,
           nonLinearGeometry::nonLinearNames_[nonLinear],
           incremental
       );

       // set traction to zero on faces not in contact
       // and set tensile stresses to zero
       //const scalar maxMagSlavePressure = gMax(mag(slavePressure_));
       scalar maxMagSlavePressure = 0.0;
       if (slavePressure_.size() > 0)
       {
           maxMagSlavePressure = max(mag(slavePressure_));
       }
       reduce(maxMagSlavePressure, maxOp<scalar>());

       forAll(touchFraction_, facei)
       {
           if
           (
               touchFraction_[facei] < SMALL
            || (
                  (slaveFaceNormals[facei] & slavePressure_[facei])
                > 1e-3*maxMagSlavePressure
            )
           )
           {
               slavePressure_[facei] = vector::zero;
               slaveValueFrac_[facei] = symmTensor::zero;
               slaveDisp_[facei] = slaveFaceNormals[facei]
                   *(slaveFaceNormals[facei]&oldSlaveDisp[facei]);
           }
       }

       // relax traction
       slavePressure_ =
           relaxFactor_*slavePressure_ + (1-relaxFactor_)*oldSlavePressure_;

       // remove any shears
       slavePressure_ = slaveFaceNormals*(slaveFaceNormals & slavePressure_);

       // limit pressure to help convergence
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

       // update old slave traction
       oldSlavePressure_ = slavePressure_;
   }

   // in parallel, the log is poluted with warnings that
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

   // under-relax value fraction
   slaveValueFrac_ =
       relaxFactor_*slaveValueFrac_ + (1-relaxFactor_)*oldSlaveValueFrac_;
   oldSlaveValueFrac_ = slaveValueFrac_;

   //--------MASTER PROCS WRITES CONTACT INFO FILE----------//
   if (Pstream::master() && (contactIterNum_ %  infoFreq_ == 0))
   {
       OFstream& contactFile = *contactFilePtr_;
       int width = 20;
       contactFile << mesh.time().value();
       contactFile.width(width);
       contactFile << contactIterNum_;
       contactFile.width(width);
       contactFile << numSlaveContactFaces;
       contactFile.width(width);
       contactFile << numSlaveSettledFaces;
       contactFile.width(width);
       contactFile << minSlavePointPenetration;
       contactFile.width(width);
       contactFile << maxMagMasterTraction;
       contactFile.width(width);
       contactFile << numCorrectedPoints;
       if (aitkenRelaxation_)
       {
           contactFile.width(width);
           contactFile << aitkenTheta_;
       }
       contactFile << endl;
   }
  //Info << "\tdone" << endl;

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


  void dirichletNeumann::writeDict(Ostream& os) const
  {
    word keyword(name()+"NormalModelDict");
    os.writeKeyword(keyword)
      << normalContactModelDict_;
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
