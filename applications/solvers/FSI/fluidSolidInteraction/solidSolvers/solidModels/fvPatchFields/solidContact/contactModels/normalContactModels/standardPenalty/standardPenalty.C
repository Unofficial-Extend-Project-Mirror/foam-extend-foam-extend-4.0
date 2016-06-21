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
    standardPenalty

\*---------------------------------------------------------------------------*/

#include "standardPenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "constitutiveModel.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(standardPenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, standardPenalty, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

void standardPenalty::calcSlavePointPoints(const label slavePatchID)
{
    // Calculate patch pointPoints
    // code adapted from enrichedPatch calcPointPoints

    // Calculate point-point addressing
    if (slavePointPointsPtr_)
    {
        FatalErrorIn
        (
            "void contactPatchPair::calcSlavePointPoints() const"
        )
            << "Point-point addressing already calculated"
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


void standardPenalty::calcPenaltyFactor()
{
    // set penalty factor
    // approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale
    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();
    const constitutiveModel& mechanical =
        mesh_.objectRegistry::lookupObject<constitutiveModel>
        (
            "rheologyProperties"
        );
    scalarField masterMu =
        mechanical.mu()().boundaryField()[masterPatchIndex];
    scalarField slaveMu =
        mechanical.mu()().boundaryField()[slavePatchIndex];
    scalarField masterLambda =
        mechanical.lambda()().boundaryField()[masterPatchIndex];
    scalarField slaveLambda =
        mechanical.lambda()().boundaryField()[slavePatchIndex];

    // Avarage contact patch bulk modulus
    scalar masterK = gAverage(masterLambda + (2.0/3.0)*masterMu);
    scalar slaveK = gAverage(slaveLambda + (2.0/3.0)*slaveMu);
    scalar bulkModulus = 0.5*(masterK+slaveK);

    // Average contact patch face area
    scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    scalar faceArea = 0.5*(masterMagSf+slaveMagSf);

    // Average contact patch cell volume
    scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    const volScalarField::DimensionedInternalField& V = mesh_.V();
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
    scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

    // Approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    penaltyFactorPtr_ =
        new scalar(penaltyScale_*bulkModulus*faceArea/cellVolume);

    Info<< "    normal penalty factor: " << *penaltyFactorPtr_ << endl;
}


void standardPenalty::calcPenetrationLimit() const
{
    if (penetrationLimitPtr_)
    {
        FatalErrorIn("void standardPenalty::calcPenetrationLimit()")
            << "penetrationLimit pointer already set" << abort(FatalError);
    }

    // We will set the penetration limit to be twice the average width of a face
    // on the contact patch

    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();

    // Average contact patch cell volume
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

    // Average cell volume
    const scalar avCellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

    // Average cell width
    const scalar avCellWidth = Foam::cbrt(avCellVolume);

    penetrationLimitPtr_ = new scalar(-2.0*avCellWidth);
}


scalar standardPenalty::penetrationLimit() const
{
    if (!penetrationLimitPtr_)
    {
        calcPenetrationLimit();
    }

    return *penetrationLimitPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPenalty::standardPenalty
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
    // reading slave pressure on restart doesn't seem to affect convergence so
    // we won't
    // (
    //     normalContactModelDict_.found("slavePressure")
    //     ? vectorField(normalContactModelDict_.lookup("slavePressure"))
    //    : vectorField(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero)
    // ),
    slaveValueFrac_
    (
        mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero
    ),
    globalSlavePointPenetration_(slaveFaceZonePatch.points().size(), 0.0),
    slavePointPenetration_
    (
        mesh_.boundaryMesh()[slavePatchID].localPoints().size(), 0.0
    ),
    areaInContact_(slaveDisp_.size(), 0.0),
    limitPenetration_
    (
        normalContactModelDict_.lookupOrDefault<Switch>
        (
            "limitPenetration", Switch(true) // default to true
        )
    ),
    penetrationLimitPtr_(NULL),
    distanceMethod_
    (
        normalContactModelDict_.lookupOrDefault<word>
        (
            "distanceMethod", "pointGGI"
        )
    ),
    correctMissedVertices_
    (
        normalContactModelDict_.lookupOrDefault<Switch>
        (
            "correctMissedVertices", Switch(false)
        )
    ),
    slavePointPointsPtr_(NULL),
    penaltyFactorPtr_(NULL),
    penaltyScale_(readScalar(normalContactModelDict_.lookup("penaltyScale"))),
    relaxFac_(readScalar(normalContactModelDict_.lookup("relaxationFactor"))),
    totalSlavePointTrac_(slavePointPenetration_.size(), 0.0),
    contactIterNum_(0),
    oscillationCorr_(normalContactModelDict_.lookup("oscillationCorrection")),
    smoothingSteps_
    (
        normalContactModelDict_.lookupOrDefault<int>("smoothingSteps", 1)
    ),
    infoFreq_(normalContactModelDict_.lookupOrDefault<int>("infoFrequency", 1)),
    contactFilePtr_(NULL)
{
    if (Pstream::master())
    {
        word masterName = mesh_.boundary()[masterPatchID].name();
        word slaveName = mesh_.boundary()[slavePatchID].name();
        fileName contactFileDir = "contact";
        mkDir(contactFileDir);
        contactFilePtr_ =
            new OFstream
            (
                contactFileDir/"normalContact_"+masterName+"_"+slaveName+".txt"
            );
        OFstream& contactFile = *contactFilePtr_;
        int width = 20;
        contactFile
            << "Time";
        contactFile.width(width);
        contactFile
        << "ContactIter";
        contactFile.width(width);
        contactFile
            << "slaveContactVer";
        contactFile.width(width);
        contactFile
            << "penetration";
        contactFile.width(width);
        contactFile
            << "maxSlaveTrac";
        contactFile.width(width);
        contactFile
            << "numMissedPoints" << endl;
    }

    // Calc point points for missed vertices
    if (correctMissedVertices_)
    {
        calcSlavePointPoints(slavePatchID);
    }

    // Calc penetration limit
    // Force calculation here as parallel operations are required
    if (limitPenetration_)
    {
        penetrationLimit();

        Info<< "        penetration limit is " << penetrationLimit() << endl;
    }
}


standardPenalty::standardPenalty(const standardPenalty& nm)
:
    normalContactModel(nm),
    normalContactModelDict_(nm.normalContactModelDict_),
    mesh_(nm.mesh_),
    slaveDisp_(nm.slaveDisp_),
    slavePressure_(nm.slavePressure_),
    slaveValueFrac_(nm.slaveValueFrac_),
    globalSlavePointPenetration_(nm.globalSlavePointPenetration_),
    slavePointPenetration_(nm.slavePointPenetration_),
    areaInContact_(nm.areaInContact_),
    limitPenetration_(nm.limitPenetration_),
    penetrationLimitPtr_(NULL),
    distanceMethod_(nm.distanceMethod_),
    correctMissedVertices_(nm.correctMissedVertices_),
    slavePointPointsPtr_(NULL),
    penaltyFactorPtr_(NULL),
    penaltyScale_(nm.penaltyScale_),
    relaxFac_(nm.relaxFac_),
    totalSlavePointTrac_(nm.totalSlavePointTrac_),
    contactIterNum_(nm.contactIterNum_),
    oscillationCorr_(nm.oscillationCorr_),
    smoothingSteps_(nm.smoothingSteps_),
    infoFreq_(nm.infoFreq_),
    contactFilePtr_(NULL)
{
    if (nm.penetrationLimitPtr_)
    {
        penetrationLimitPtr_ = new scalar(*nm.penetrationLimitPtr_);
    }

    if (nm.slavePointPointsPtr_)
    {
        slavePointPointsPtr_ = new labelListList(*nm.slavePointPointsPtr_);
    }

    if (nm.penaltyFactorPtr_)
    {
        penaltyFactorPtr_ = new scalar(*nm.penaltyFactorPtr_);
    }

    if (nm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*nm.contactFilePtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //


void standardPenalty::correct
(
    const primitiveFacePatch& masterFaceZonePatch,
    const primitiveFacePatch& slaveFaceZonePatch,
    const intersection::algorithm alg,
    const intersection::direction dir,
    word fieldName,
    const Switch orthotropic,
    const nonLinearGeometry::nonLinearType nonLinear,
    vectorField& slaveFaceNormals,
    const ggiZoneInterpolation& zoneToZone
)
{
    //- Preliminaries
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    scalar maxMagSlaveTraction = 0.0;
    contactIterNum_++;

    //--------------------CALCULATE SLAVE PENETRATION----------------------//

    if (distanceMethod_ == "pointGGI")
    {
        // New distance method 'pointGGI" is fastest
        globalSlavePointPenetration_ =
            zoneToZone.slavePointDistanceToIntersection();
    }
    else if (distanceMethod_ == "point")
    {
        // Create master to slave interpolation for calculation of distances
        PatchToPatchInterpolation<primitiveFacePatch, primitiveFacePatch>
            masterToSlavePatchToPatchInterpolator
            (
                masterFaceZonePatch, // from zone
                slaveFaceZonePatch, // to zone
                alg,
                dir
            );

        // patchToPatch
        globalSlavePointPenetration_ =
            masterToSlavePatchToPatchInterpolator.pointDistanceToIntersection();
    }
    else
    {
        FatalError
            << "distanceMethod " << distanceMethod_ << " not allowed" << nl
            << "valid options are: pointGGI and point" << abort(FatalError);
    }

    forAll(globalSlavePointPenetration_, pointI)
    {
        // When the master surface surrounds the slave (like the pelvis and
        // femur head) then the slave penetration can sometimes calculate the
        // distance through the femur head to the pelvis which is wrong so I
        // limit slavePenetration here
        if (limitPenetration_)
        {
            if (globalSlavePointPenetration_[pointI] < penetrationLimit())
            {
                globalSlavePointPenetration_[pointI] = 0.0;
            }
        }
    }

    forAll(slavePointPenetration_, pointI)
    {
        // The local point values are kept at the start of the global field
        slavePointPenetration_[pointI] =
            globalSlavePointPenetration_
            [
                pointI
            ];
    }

    label numCorrectedPoints = 0;
    if (correctMissedVertices_)
    {
#       include "iterativePenaltyCorrectMissedVertices.H"
    }

    // Calculate slaveFaceNormals
    // We have a number of options:
    //     1) use the slave deformed/undeformed normals;
    //     2) use the master deformed/undeformed normals interpolated to the
    //        slave;
    //     3) use deformed normals calculated using deformation gradient
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
    forAll(slavePatchLocalFaces, facei)
    {
        areaInContact_[facei] =
            slavePatchLocalFaces[facei].areaInContact
            (
                slavePatchLocalPoints,
                slavePointPenetration_
            );
    }

    int numSlaveContactPoints = 0;

    scalar minSlavePointPenetration = gMin(globalSlavePointPenetration_);

    scalar penaltyFac = penaltyFactor();

    // Calculate traction increments
    forAll(totalSlavePointTrac_, pointI)
    {
        // If a point has penetrated (i.e. if the penetration is negative),
        if (slavePointPenetration_[pointI] < 0.0)
        {
            // Count points in contact
            numSlaveContactPoints++;

            // The force is linearly proportional the penetration, like a spring
            totalSlavePointTrac_[pointI] =
                penaltyFac*slavePointPenetration_[pointI];
        }
        else
        {
            totalSlavePointTrac_[pointI] = 0.0;
        }
    }

    //--------INTERPOLATE SLAVE POINT TRACTION TO SLAVE FACES----------//

    // Create local patch interpolation: No need to interpolate using the entire
    // face zone patch
    primitivePatchInterpolation localSlaveInterpolator
        (
            mesh.boundaryMesh()[slavePatchIndex]
        );

    vectorField newSlavePressure =
        localSlaveInterpolator.pointToFaceInterpolate<scalar>
        (
            totalSlavePointTrac_
        )*slaveFaceNormals;

    // Oscillations sometimes appear in the contact pressure
    // We will perform some simple smoothing here to quell them
    if (oscillationCorr_)
    {
        if (smoothingSteps_ < 1)
        {
            FatalError
                << "smoothingSteps must be greater than or equal to 1 if"
                << " oscillation correction is active" << abort(FatalError);
        }

        // Interpolate face values to points then interpolate back which
        // smoothes the field
        primitivePatchInterpolation localSlaveInterpolator
            (
                mesh.boundaryMesh()[slavePatchIndex]
            );

        vectorField newSlavePressurePoints
            (
                mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero
            );

        // Contact area tends to increase with smoothing so we will mark the
        // faces in contact before smoothing, and eliminate traction on all
        // other faces after smoothing
        boolList facesInContact(newSlavePressure.size(), false);
        forAll(facesInContact, faceI)
        {
            if (mag(newSlavePressure[faceI]) > SMALL)
            {
                facesInContact[faceI] = true;
            }
        }

        for (int i=0; i<smoothingSteps_; i++)
        {
            newSlavePressurePoints =
                localSlaveInterpolator.faceToPointInterpolate<vector>
                (
                    newSlavePressure
                );
            newSlavePressure =
                localSlaveInterpolator.pointToFaceInterpolate<vector>
                (
                    newSlavePressurePoints
                );
        }

        // Remove faces introduced to the contact are due to smoothing
        forAll(facesInContact, faceI)
        {
            if (!facesInContact[faceI])
            {
                newSlavePressure[faceI] = vector::zero;
            }
        }
    }

    // Under-relax pressure
    slavePressure_ =
        relaxFac_*newSlavePressure + (1.0 - relaxFac_)*slavePressure_;

    // Add slavePressure to dict to allow restart
    // PC 19-Nov-14: disabled
    //normalContactModelDict_.set("slavePressure", slavePressure_);

    // In parallel, the log is poluted with warnings about getting the max of a
    // list of size zero so we will only use procs with slave faces when
    // calculating the max
    // maxMagSlaveTraction = gMax(mag(slavePressure_));
    if (slavePressure_.size() > 0)
    {
        maxMagSlaveTraction = max(mag(slavePressure_));
    }
    reduce(maxMagSlaveTraction, maxOp<scalar>());

    reduce(numSlaveContactPoints, sumOp<int>());

    // Write to contact file
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
            << numSlaveContactPoints;
        contactFile.width(width);
        contactFile
            << minSlavePointPenetration;
        contactFile.width(width);
        contactFile
            << maxMagSlaveTraction;
        contactFile.width(width);
        contactFile
            << numCorrectedPoints << endl;
    }
}


void standardPenalty::resetSlavePressure()
{
    slavePressure_ = vector::zero;

    normalContactModelDict_.remove("slavePressure");
}


void standardPenalty::writeDict(Ostream& os) const
{
    word keyword(name() + "NormalModelDict");
    os.writeKeyword(keyword)
        << normalContactModelDict_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
