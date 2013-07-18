/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

Description
    A pair of surfaces in contact.

\*---------------------------------------------------------------------------*/

#include "contactPatchPair.H"
#include "contactProblem.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::contactPatchPair::contactPatchPair
(
    const word& name,
    contactProblem& cp,
    const dictionary& dict
)
:
    contactActive_(dict.lookup("contactActive")),
    name_(name),
    cp_(cp),
    alg_(intersection::algorithmNames_.read(dict.lookup("projectionAlgo"))),
    dir_(intersection::directionNames_.read(dict.lookup("projectionDir"))),
    masterPatch_(dict.lookup("masterPatch"), cp.mesh().boundaryMesh()),
    slavePatch_(dict.lookup("slavePatch"), cp.mesh().boundaryMesh()),
    contactGapTol_(readScalar(dict.lookup("contactGapTol"))),
    penaltyScaleFactor_(readScalar(dict.lookup("penaltyScaleFactor"))),
    returnScaleFactor_(readScalar(dict.lookup("returnScaleFactor"))),
    contactIterNum_(0),
    contactFilePtr_(NULL),
    slaveFaceZoneName_(slavePatch_.name() + "FaceZone"),
    masterFaceZoneName_(masterPatch_.name() + "FaceZone"),
    slavePointZoneName_(slavePatch_.name() + "PointZone"),
    totalSlavePointForce_(
  			  cp_.mesh().boundaryMesh()[slavePatch_.index()].pointNormals().size(), vector::zero),
			  // IOobject
			  // (
			  //  "totalSlavePointForce",
			  //  cp_.mesh().time().timeName(),
			  //  cp_.mesh(),
			  //  IOobject::READ_IF_PRESENT,
			  //  IOobject::AUTO_WRITE
			  //  ),
    			  // vectorField(cp_.mesh().boundaryMesh()[slavePatch_.index()].pointNormals().size(), vector::zero)
			  // ),
    slavePointPenetration_(
			   cp_.mesh().boundaryMesh()[slavePatch_.index()].pointNormals().size(),
			   0.0
			   ),
    masterPointPenetration_(
			    cp_.mesh().boundaryMesh()[masterPatch_.index()].pointNormals().size(),
			    0.0
			    ),
    globalSlavePointPenetration_(
				 cp_.mesh().pointZones()[cp_.mesh().faceZones().findZoneID(slaveFaceZoneName_)].size(),
				 0.0
				 ),
    globalMasterPointPenetration_(
				  cp_.mesh().pointZones()[cp_.mesh().faceZones().findZoneID(masterFaceZoneName_)].size(),
				  0.0
				  ),
    oldTotalSlavePointForce_(
    			     cp_.mesh().boundaryMesh()[slavePatch_.index()].pointNormals().size(),
    			     vector::zero
    			     ),
    oldMinSlavePointPenetration_(0.0),
    penetrationLimit_(readScalar(dict.lookup("penetrationLimit"))),
    rigidMaster_(dict.lookup("rigidMaster")),
    interpolationMethod_(dict.lookup("interpolationMethod")),
    faceZoneMasterToSlaveInterpolator_(
				       cp_.mesh().faceZones()[cp_.mesh().faceZones().findZoneID(masterFaceZoneName_)](), // from
				       cp_.mesh().faceZones()[cp_.mesh().faceZones().findZoneID(slaveFaceZoneName_)](), // to zone
				       alg_,
				       dir_
				       ),
    faceZoneSlaveToMasterInterpolatorPtr_(NULL),
    faceZoneGgiInterpolatorPtr_(NULL),
    slaveInterpolator_(cp_.mesh().boundaryMesh()[slavePatch_.index()])
{
  Info << "\tConstructing contact patch pair"
       << "\n\t\tmaster patch:\t" << masterPatch_.name()
       << "\n\t\tslave patch:\t" << slavePatch_.name()
       << endl;

  if(rigidMaster_)
    Info << "\t\tThe master surface is considered rigid and is set as traction-free" << endl;

  //- Check interpolation scheme for passing tractions from slave to master
  if(interpolationMethod_ == "patchToPatch")
    {
      Info << "\t\tMethod for interpolation of traction from slave to master: patchToPatch" << endl;
      label masterFaceZoneID = cp_.mesh().faceZones().findZoneID(masterFaceZoneName_);
      label slaveFaceZoneID = cp_.mesh().faceZones().findZoneID(slaveFaceZoneName_);
      faceZoneSlaveToMasterInterpolatorPtr_ =
	new zoneToZoneInterpolation
	(
	 cp_.mesh().faceZones()[slaveFaceZoneID](), // from zone
	 cp_.mesh().faceZones()[masterFaceZoneID](), // to zone
	 alg_,
	 dir_
	 );
    }
  else if(interpolationMethod_ == "ggi")
    {
      Info << "\t\tMethod for interpolation of traction from slave to master: ggi" << endl;
      label masterFaceZoneID = cp_.mesh().faceZones().findZoneID(masterFaceZoneName_);
      label slaveFaceZoneID = cp_.mesh().faceZones().findZoneID(slaveFaceZoneName_);
      faceZoneGgiInterpolatorPtr_ =
	new ggiZoneInterpolation
	(
	 cp_.mesh().faceZones()[masterFaceZoneID](), // master zone
	 cp_.mesh().faceZones()[slaveFaceZoneID](), // slave zone
	 tensorField(0),
	 tensorField(0),
	 vectorField(0),
	 0.0,
	 0.0,
	 true,
	 ggiInterpolation::AABB
	 );
    }
  else
    {
      SeriousError << "\n\nTraction interpolation method '"
		   << interpolationMethod_ << "' not found!\n"
		   << "\nValid interpolation methods are:\n"
		   << "ggi\n"
		   << "patchToPatch"
		   << endl
		   << exit(FatalError);
    }

  //- only the master should create a contactInfo file
  if(Pstream::master())
    {
      contactFilePtr_ = new OFstream(name_);
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
  Info << "\tContact patch pair constructed" << endl;
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::masterTouchFraction() const
{
  //--DOES NOT WORK IN PARALLEL YET BUT IS ONLY FOR VISUALISATION--/

  const fvMesh& mesh = cp_.mesh();

  scalarField pDistToInter(cp_.mesh().boundaryMesh()[masterPatch_.index()].size(), 0.0);
  //    = slaveToMasterInterpolateDeformed.pointDistanceToIntersection();

  scalarField vertexMasterGap =  pDistToInter;

// Calculate area in contact
    const faceList& masterPatchLocalFaces =
        mesh.boundaryMesh()[masterPatch_.index()].localFaces();

    const pointField& masterPatchLocalPoints =
        mesh.boundaryMesh()[masterPatch_.index()].localPoints();

    tmp<scalarField> ttouchFrac
    (
        new scalarField(masterPatchLocalFaces.size(), 0)
    );
    scalarField& touchFrac = ttouchFrac();

    forAll (masterPatchLocalFaces, faceI)
    {
        touchFrac[faceI] =
            masterPatchLocalFaces[faceI].areaInContact
            (
                masterPatchLocalPoints,
                vertexMasterGap
            );
    }

    return ttouchFrac;
}




Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::slaveTouchFraction() const
{
  //--DOES NOT WORK IN PARALLEL YET BUT IS ONLY FOR VISUALISATION--/

  const fvMesh& mesh = cp_.mesh();

  scalarField vertexSlaveGap = slavePointPenetration_;
  //    (cp_.mesh().boundaryMesh()[slavePatch_.index()].size(), 0.0);

  // Calculate area in contact
    const faceList& slavePatchLocalFaces =
      mesh.boundaryMesh()[slavePatch_.index()].localFaces();

    const pointField& slavePatchLocalPoints =
        mesh.boundaryMesh()[slavePatch_.index()].localPoints();

    tmp<scalarField> ttouchFrac
    (
        new scalarField(slavePatchLocalFaces.size(), 0)
    );
    scalarField& touchFrac = ttouchFrac();

    forAll (slavePatchLocalFaces, faceI)
    {
        touchFrac[faceI] =
            slavePatchLocalFaces[faceI].areaInContact
            (
                slavePatchLocalPoints,
                vertexSlaveGap
            );
    }

    return ttouchFrac;
}




Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::masterGapPoints() const
{
  scalarField globalMasterPointPenetration
    (
     cp_.mesh().boundaryMesh()[masterPatch_.index()].meshPoints().size(),
     0.0
     );

  scalarField vertexMasterGap = masterPointPenetration_;
  //globalMasterPointPenetration;

    tmp<scalarField> tcontactGapPoints
    (
	  new scalarField(vertexMasterGap.size(), 0)
    );
    scalarField& contactGapPoints = tcontactGapPoints();

    forAll (vertexMasterGap, pointI)
    {
      contactGapPoints[pointI] = vertexMasterGap[pointI];
    }

    return tcontactGapPoints;
}




Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::slaveGapPoints() const
{
  const scalarField& slavePointPenetration = slavePointPenetration_;

  tmp<scalarField> tcontactGapPoints
    (
     new scalarField(slavePointPenetration.size(), 0)
    );
    scalarField& contactGapPoints = tcontactGapPoints();

    forAll (slavePointPenetration, pointI)
    {
      contactGapPoints[pointI] = slavePointPenetration[pointI];
    }

    return tcontactGapPoints;
}




Foam::tmp<Foam::vectorField>
Foam::contactPatchPair::masterPointForce() const
{
  //- returns zero at the moment
  vectorField masterPointForce
    (
     cp_.mesh().boundaryMesh()[masterPatch_.index()].meshPoints().size(),
     vector::zero
     );

  tmp<vectorField> tcontactPointForce
    (
     new vectorField(masterPointForce.size(), vector::zero)
     );
  vectorField& contactPointForce = tcontactPointForce();

  forAll (contactPointForce, pointI)
    {
      contactPointForce[pointI] = masterPointForce[pointI];
    }

  return tcontactPointForce;
}




Foam::tmp<Foam::vectorField>
Foam::contactPatchPair::slavePointForce() const
{
  vectorField slavePointForce = totalSlavePointForce_;

  tmp<vectorField> tcontactPointForce
    (
     new vectorField(slavePointForce.size(), vector::zero)
     );

  vectorField& contactPointForce = tcontactPointForce();

  forAll (contactPointForce, pointI)
    {
      contactPointForce[pointI] = slavePointForce[pointI];
    }

  return tcontactPointForce;
}




Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::masterContactPressure() const
{
    const fvMesh& mesh = cp_.mesh();

    vectorField n = mesh.boundary()[masterPatch_.index()].nf();

    const volSymmTensorField& sigma =
        mesh.objectRegistry::lookupObject<volSymmTensorField>("sigma");

    scalarField masterPressure
      = n & (n & sigma.boundaryField()[masterPatch_.index()]);

    tmp<scalarField> tcontactPressure
    (
        new scalarField(masterPressure.size(), 0)
    );
    scalarField& contactPressure = tcontactPressure();

    forAll (masterPressure, faceI)
    {
      contactPressure[faceI] = max(-1*masterPressure[faceI],0);
    }

    return tcontactPressure;
}




Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::slaveContactPressure() const
{
    const fvMesh& mesh = cp_.mesh();

    vectorField n = mesh.boundary()[slavePatch_.index()].nf();

    const volSymmTensorField& sigma =
        mesh.objectRegistry::lookupObject<volSymmTensorField>("sigma");

    scalarField slavePressure
      = n & (n & sigma.boundaryField()[slavePatch_.index()]);

    tmp<scalarField> tcontactPressure
    (
        new scalarField(slavePressure.size(), 0)
    );
    scalarField& contactPressure = tcontactPressure();

    forAll (slavePressure, faceI)
    {
      contactPressure[faceI] = max(-1*slavePressure[faceI],0);
    }

    return tcontactPressure;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
