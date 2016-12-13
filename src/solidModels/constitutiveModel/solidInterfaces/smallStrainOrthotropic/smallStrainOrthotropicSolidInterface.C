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

#include "smallStrainOrthotropicSolidInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smallStrainOrthotropicSolidInterface, 0);
    addToRunTimeSelectionTable
    (solidInterface, smallStrainOrthotropicSolidInterface, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::smallStrainOrthotropicSolidInterface::smallStrainOrthotropicSolidInterface
(
 const word& name,
 const fvMesh& mesh,
 const constitutiveModel& rheology
)
:
  solidInterface(name, mesh, rheology)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smallStrainOrthotropicSolidInterface::
~smallStrainOrthotropicSolidInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smallStrainOrthotropicSolidInterface::correct(fvVectorMatrix& UEqn)
{
  const fvMesh& mesh = solidInterface::mesh();

  const unallocLabelList& owner = mesh.owner();
  const unallocLabelList& neighbour = mesh.neighbour();

  const volVectorField& U = UEqn.psi();
  const vectorField& UI = U.internalField();

  const volTensorField& gradU =
    mesh.lookupObject<volTensorField>("grad(" + U.name() + ')');
  const tensorField& gradUI = gradU.internalField();

  const volSymmTensor4thOrderField& matC =
      mesh.lookupObject<volSymmTensor4thOrderField>("C");
  const volDiagTensorField& K = mesh.lookupObject<volDiagTensorField>("K");

  const vectorField& SI  = mesh.Sf().internalField();
  const scalarField& magSI  = mesh.magSf().internalField();
  const scalarField& deltaCoeffs = mesh.deltaCoeffs().internalField();
  const scalarField& w = mesh.weights().internalField();

  scalarField& diag = UEqn.diag();
  scalarField& upper = UEqn.upper();
  vectorField& source = UEqn.source();
  FieldField<Field, vector>& boundaryCoeffs = UEqn.boundaryCoeffs();

  vectorField& interU = interfaceDisplacement();


  // Internal faces
  forAll(globalInterFaces(), faceI)
    {
      label curGlobalFace = globalInterFaces()[faceI];

      label curOwner = owner[curGlobalFace];
      label curNeighbour = neighbour[curGlobalFace];

      vector ownN = SI[curGlobalFace]/magSI[curGlobalFace];
      //vector ngbN = -ownN;

      scalar magS = magSI[curGlobalFace];

      // the interface tangential gradient may be calculated
      // in one of three ways:
      //    1) extrapolation from adjoining cell centres
      //    2) inteprolation from adjoining cell centres
      //    3) directly calculated using the finite area method

      // extrapolate tangential gradient
      tensor ownGradU = gradUI[curOwner];
      tensor ngbGradU = gradUI[curNeighbour];

      scalar ngbDn = w[curGlobalFace]*(1.0/deltaCoeffs[curGlobalFace]);
      scalar ownDn = (1.0/deltaCoeffs[curGlobalFace]) - ngbDn;

      symmTensor4thOrder ownMatC = matC.internalField()[curOwner];
      symmTensor4thOrder ngbMatC = matC.internalField()[curNeighbour];

      diagTensor ownK = K.internalField()[curOwner];
      diagTensor ngbK = K.internalField()[curNeighbour];

      vector ownU = UI[curOwner];
      vector ngbU = UI[curNeighbour];

      // Interface displacement
      // no decomposition of traction - philipc

      scalar ownKn = ownN & (ownN & ownK);
      scalar ngbKn = ownN & (ownN & ngbK);
      vector ownKt = (I - sqr(ownN)) & (ownN & ownK);
      vector ngbKt = (I - sqr(ownN)) & (ownN & ngbK);

      symmTensor ownEpsilon = symm(ownGradU);
      symmTensor ngbEpsilon = symm(ngbGradU);

      // explicit terms
      vector ownQ =
          (ownKt & ownGradU)
          + ( ownN & ((ownMatC && ownEpsilon) - (ownK & ownGradU)) );
      vector ngbQ =
          (ngbKt & ngbGradU)
          + ( ownN & ((ngbMatC && ngbEpsilon) - (ngbK & ngbGradU)) );

      // no decomposition
      interU[faceI] =
    (
     ownKn*ngbDn*ownU
     + ngbKn*ownDn*ngbU
     + ownDn*ngbDn*(ngbQ - ownQ)
     )
    /(ownKn*ngbDn + ngbKn*ownDn);


      // Implicit coupling

      scalar wRevLin = 1.0 - w[curGlobalFace];

      scalar Knf = 1.0/(wRevLin/ownKn + (1.0-wRevLin)/ngbKn);

      scalar Dnf = 1.0/deltaCoeffs[curGlobalFace];

      // Owner
      //diag[curOwner] += Kf*magS/Dnf;
      diag[curOwner] += Knf*magS/Dnf; //Kf*magS/Dnf;

      //upper[curGlobalFace] -= Kf*magS/Dnf;
      upper[curGlobalFace] -= Knf*magS/Dnf; //Kf*magS/Dnf;

      source[curOwner] +=
    (
     ownKn*ngbDn*ngbQ
     + ngbKn*ownDn*ownQ
     )*magS
    /(ownKn*ngbDn + ngbKn*ownDn);


      // Neighbour
      diag[curNeighbour] += Knf*magS/Dnf; //Kf*magS/Dnf;

      // philipc - no decomposition
        source[curNeighbour] -=
      (
       ownKn*ngbDn*ngbQ
       + ngbKn*ownDn*ownQ
       )*magS
      /(ownKn*ngbDn + ngbKn*ownDn);
    }


  // Processor faces
  forAll(processorPatchFaces(), patchI)
    {
      label curPatch = processorPatches()[patchI];

      vectorField& curProcInterU =
    processorInterfaceDisplacement()[patchI];

      const vectorField curProcOwnU =
    U.boundaryField()[curPatch].patchInternalField();
      const vectorField curProcNgbU =
    U.boundaryField()[curPatch].patchNeighbourField();

      const tensorField curProcOwnGradU =
    gradU.boundaryField()[curPatch].patchInternalField();
      const tensorField curProcNgbGradU =
    gradU.boundaryField()[curPatch].patchNeighbourField();

      const vectorField& curProcS =
    mesh.Sf().boundaryField()[curPatch];
      const scalarField& curProcMagS =
    mesh.magSf().boundaryField()[curPatch];
      const scalarField& curProcDeltaCoeffs =
    mesh.deltaCoeffs().boundaryField()[curPatch];
      const scalarField& curProcW =
    mesh.weights().boundaryField()[curPatch];
      const vectorField curProcOwnC =
    mesh.C().boundaryField()[curPatch].patchInternalField();
      const vectorField curProcNgbC =
    mesh.C().boundaryField()[curPatch].patchNeighbourField();

      const symmTensor4thOrderField curProcOwnMatC =
    matC.boundaryField()[curPatch].patchInternalField();
      const symmTensor4thOrderField curProcNgbMatC =
    matC.boundaryField()[curPatch].patchNeighbourField();

      const diagTensorField curProcOwnK =
    K.boundaryField()[curPatch].patchInternalField();
      const diagTensorField curProcNgbK =
    K.boundaryField()[curPatch].patchNeighbourField();

      const unallocLabelList& curProcFaceCells =
    mesh.boundary()[curPatch].faceCells();

      forAll(processorPatchFaces()[patchI], faceI)
        {
      label curFace = processorPatchFaces()[patchI][faceI];

      scalar ngbDn = curProcW[curFace]*(1.0/curProcDeltaCoeffs[curFace]);
      scalar ownDn = (1.0/curProcDeltaCoeffs[curFace]) - ngbDn;

      scalar magS = curProcMagS[curFace];

      symmTensor4thOrder ownMatC  = curProcOwnMatC[curFace];
      symmTensor4thOrder ngbMatC  = curProcNgbMatC[curFace];

      diagTensor ownK  = curProcOwnK[curFace];
      diagTensor ngbK  = curProcNgbK[curFace];

      vector ownN = curProcS[curFace]/curProcMagS[curFace];
      //vector ngbN = -ownN;

      // extrapolate tangential gradient
      tensor ownGradU = curProcOwnGradU[curFace];
      tensor ngbGradU = curProcNgbGradU[curFace];

      vector ownU = curProcOwnU[curFace];
      vector ngbU = curProcNgbU[curFace];

      // Interface displacement

      // stiffnesses
      scalar ownKn = ownN & (ownN & ownK);
      scalar ngbKn = ownN & (ownN & ngbK);
      vector ownKt = (I - sqr(ownN)) & (ownN & ownK);
      vector ngbKt = (I - sqr(ownN)) & (ownN & ngbK);

      symmTensor ownEpsilon = symm(curProcOwnGradU[curFace]);
      symmTensor ngbEpsilon = symm(curProcNgbGradU[curFace]);

      // explicit terms
      vector ownQ =
          (ownKt & ownGradU)
          + ( ownN & ((ownMatC && ownEpsilon) - (ownK & ownGradU)) );
      vector ngbQ =
          (ngbKt & ngbGradU)
          + ( ownN & ((ngbMatC && ngbEpsilon) - (ngbK & ngbGradU)) );

      curProcInterU[faceI] =
        (
         ownKn*ngbDn*ownU
         + ngbKn*ownDn*ngbU
         + ownDn*ngbDn*(ngbQ - ownQ)
         )
        /(ownKn*ngbDn + ngbKn*ownDn);


      // Implicit coupling

      scalar wRevLin = 1.0 - curProcW[curFace];

      scalar Knf = 1.0/(wRevLin/ownKn + (1.0-wRevLin)/ngbKn);

      scalar Dnf = 1.0/curProcDeltaCoeffs[curFace];

      // Owner
      diag[curProcFaceCells[curFace]] += Knf*magS/Dnf;

      boundaryCoeffs[curPatch][curFace] += Knf*magS*vector::one/Dnf;

      // no decomposition
      source[curProcFaceCells[curFace]] +=
        (
         ownKn*ngbDn*ngbQ
         + ngbKn*ownDn*ownQ
         )*magS
        /(ownKn*ngbDn + ngbKn*ownDn);
        }
    }
}

// ************************************************************************* //
