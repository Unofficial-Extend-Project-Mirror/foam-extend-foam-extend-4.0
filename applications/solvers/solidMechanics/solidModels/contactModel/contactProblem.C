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
    Class describes a multiple body contact problem.  Each individual contact
    is described by a contactPatchPair.  contactProblem handles
    multiple contact updates and sets the boundary conditions on the
    displacement field.

\*---------------------------------------------------------------------------*/

#include "contactProblem.H"
#include "fvMesh.H"
#include "FieldFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "surfaceFields.H"
#include "pointMesh.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactProblem, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Read constructor given IOobject
contactProblem::contactProblem
(
 volVectorField& U
)
:
    IOdictionary
    (
        IOobject
        (
            "contactProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    contactPatchPairList(),
    U_(U)
{
  Info << "\nConstructing contact problem" << endl;
  Info << "\t*************************************************************************************\n"
       << "\t** MAKE SURE MASTER AND SLAVE FACE AND POINT ZONES HAVE BEEN DEFINED               **\n"
       << "\t** To define, use the 'setSet' utility:                                            **\n"
       << "\t** faceSet <slaveName>FaceZone new patchToFace <slaveName>                         **\n"
       << "\t** faceSet <masterName>FaceZone new patchToFace <masterName>                       **\n"
       << "\t** pointSet <masterName>PointZone new faceToPoint <masterName>FaceZone all         **\n"
       << "\t** pointSet <slaveName>PointZone new faceToPoint <slaveName>FaceZone all           **\n"
       << "\t** Then use the 'setsToZone -noFlipMap' command                                    **\n"
       << "\t** For parallel runs, 'globalFaceZones (<slaveName>FaceZone <masterName>FaceZone)' **\n"
       << "\t** must be included in the decomposeParDict                                        **\n"
       << "\t** <slaveName> and <masterName> are replaced with the slave and master patch names **\n"
       << "\t*************************************************************************************"
       << endl;

  //- Read contactPatchPairList
  Istream& is = lookup("contacts");

  PtrList<entry> contactEntries(is);
  
  contactPatchPairList& contacts = *this;
  
  contacts.setSize(contactEntries.size());
  
  forAll(contacts, contactI)
    {
      contacts.set
        (
	 contactI,
	 new contactPatchPair
	 (
	  contactEntries[contactI].keyword(),
	  *this,
	  contactEntries[contactI].dict()
	  )
	 );
     }

  Info << "Contact problem constructed"
       << endl;
}
  

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//**********************CORRECT FUNCTION*************************************//
void contactProblem::correct()
{
  contactPatchPairList& contacts = *this;

  // Collect patches involved in contact
  boolList contactPatches(U().boundaryField().size(), false);

  forAll (contacts, contactI)
    {
      contactPatches[contacts[contactI].masterPatch().index()] = true;
      contactPatches[contacts[contactI].slavePatch().index()] = true;
    }

  // Calculate contact trcations
  forAll (contacts, contactI)
    {
      if(contacts[contactI].contactActive())
	{
	  contacts[contactI].correct();
	}
      else
	{
	  Info << "\t\t\tContact " << contacts[contactI].name() << " not active" << endl;
	}
    }
}




//**********************CONTACT AREA FUNCTION***********************************//
tmp<volScalarField> contactProblem::contactArea() const
{
    tmp<volScalarField> tca
    (
        new volScalarField
        (
            IOobject
            (
                "contactArea",
                U().time().timeName(),
                U().db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(0)
        )
    );

    volScalarField& ca = tca();

    // Set contact area boundary
    const contactPatchPairList& contacts = *this;

    forAll (contacts, contactI)
    {
        // Get master contact
        ca.boundaryField()[contacts[contactI].masterPatch().index()] +=
            contacts[contactI].masterTouchFraction();

        // Get slave contact
        ca.boundaryField()[contacts[contactI].slavePatch().index()] +=
            contacts[contactI].slaveTouchFraction();

    //-----------CALCULATE-ACTUAL-CONTACT-AREA--------------//
    label masterIndex = contacts[contactI].masterPatch().index();
    label slaveIndex = contacts[contactI].slavePatch().index();
    scalarField masterFrac = contacts[contactI].masterTouchFraction();
    scalarField slaveFrac = contacts[contactI].slaveTouchFraction();
    scalar contactAreaMaster =
      gSum
      (
       masterFrac *
       mag(
	   mesh().Sf().boundaryField()[masterIndex]
	   )
       );
    scalar contactAreaSlave =
      gSum
      (
       slaveFrac *
       mag(
	   mesh().Sf().boundaryField()[slaveIndex]
	   )
       );  
    Info << "\nContact area of master patch is: "
	 << contactAreaMaster << " m^2"
	 << "\nContact area of slave patch is: "
	 << contactAreaSlave << " m^2"
	 << endl << endl;
    //------------------------------------------------------//
    }

    return tca;
}




//tmp<pointScalarField> contactProblem::contactGapPoints() const
void contactProblem::contactGapPoints(pointScalarField& cGapPoints)
{
  const  contactPatchPairList& contacts = *this;
  
  scalarField& cGapPointsInternal = cGapPoints.internalField();
  
  forAll (contacts, contactI)
    {
      scalarField masterGapPoints = contacts[contactI].masterGapPoints();
      labelList masterBoundaryLabels = mesh().boundaryMesh()[contacts[contactI].masterPatch().index()].meshPoints();       
      
      scalarField slaveGapPoints = contacts[contactI].slaveGapPoints();
      labelList slaveBoundaryLabels = mesh().boundaryMesh()[contacts[contactI].slavePatch().index()].meshPoints();       

      forAll(masterBoundaryLabels, pointI)
	{
	  cGapPointsInternal[masterBoundaryLabels[pointI]] = masterGapPoints[pointI];
	}
      forAll(slaveBoundaryLabels, pointI)
	{
	  cGapPointsInternal[slaveBoundaryLabels[pointI]] = slaveGapPoints[pointI];
	}
    }
}




void contactProblem::contactPointForce(pointVectorField& cPointForce)
{
  pointMesh pMesh(mesh());
  const contactPatchPairList& contacts = *this;
  
  vectorField& cPointForceInternal = cPointForce.internalField();
  
  forAll (contacts, contactI)
    {
      vectorField masterContactPointForce = contacts[contactI].masterPointForce();
      labelList masterBoundaryLabels = pMesh.boundary()[contacts[contactI].masterPatch().index()].meshPoints();
      
      vectorField slaveContactPointForce = contacts[contactI].slavePointForce();
      labelList slaveBoundaryLabels = pMesh.boundary()[contacts[contactI].slavePatch().index()].meshPoints();
      
      forAll(masterBoundaryLabels, pointI)
	{
	  cPointForceInternal[masterBoundaryLabels[pointI]] = masterContactPointForce[pointI];
	}
      forAll(slaveBoundaryLabels, pointI)
	{
	  cPointForceInternal[slaveBoundaryLabels[pointI]] = slaveContactPointForce[pointI];
	}
    }
}




tmp<volScalarField> contactProblem::contactPressure() const
{
    tmp<volScalarField> tcPress
    (
        new volScalarField
        (
            IOobject
            (
                "contactPressure",
                U().time().timeName(),
                U().db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(0)
        )
    );

    volScalarField& cPress = tcPress();

    const contactPatchPairList& contacts = *this;

    forAll (contacts, contactI)
    {
        // Get master contact pressure
        cPress.boundaryField()[contacts[contactI].masterPatch().index()] +=
            contacts[contactI].masterContactPressure();

        // Get slave contact pressure
        cPress.boundaryField()[contacts[contactI].slavePatch().index()] +=
            contacts[contactI].slaveContactPressure();
    }

    return tcPress;
}

// Return a list of contactPatchPair names
wordList contactProblem::names() const
{
    const contactPatchPairList& contacts = *this;

    wordList t(contacts.size());

    forAll (contacts, contactI)
    {
        t[contactI] = contacts[contactI].name();
    }

    return t;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
