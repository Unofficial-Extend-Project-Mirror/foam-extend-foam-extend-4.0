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

#include "smallStrainCorrectedSolidInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smallStrainCorrectedSolidInterface, 0);
    addToRunTimeSelectionTable
    (solidInterface, smallStrainCorrectedSolidInterface, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::smallStrainCorrectedSolidInterface::smallStrainCorrectedSolidInterface
(
 const word& name,
 const fvMesh& mesh,
 const constitutiveModel& rheology
)
:
  solidInterface(name, mesh, rheology)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smallStrainCorrectedSolidInterface::~smallStrainCorrectedSolidInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smallStrainCorrectedSolidInterface::correct(fvVectorMatrix& UEqn)
{
  const fvMesh& mesh = solidInterface::mesh();

  const unallocLabelList& owner = mesh.owner();
  const unallocLabelList& neighbour = mesh.neighbour();

  const volVectorField& U = UEqn.psi();
  const vectorField& UI = U.internalField();

  const volTensorField& gradU =
    mesh.lookupObject<volTensorField>("grad(" + U.name() + ')');
  const tensorField& gradUI = gradU.internalField();

  const volScalarField& mu = mesh.lookupObject<volScalarField>("mu");
  const volScalarField& lambda = mesh.lookupObject<volScalarField>("lambda");

  const vectorField& SI  = mesh.Sf().internalField();
  const scalarField& magSI  = mesh.magSf().internalField();
  const scalarField& deltaCoeffs = mesh.deltaCoeffs().internalField();
  const scalarField& w = mesh.weights().internalField();
  const vectorField& CI  = mesh.C().internalField();
  const vectorField& CfI  = mesh.Cf().internalField();

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

      scalar ownTrGradU = tr(ownGradU);
      scalar ngbTrGradU = tr(ngbGradU);

      scalar ngbDn = w[curGlobalFace]*(1.0/deltaCoeffs[curGlobalFace]);
      scalar ownDn = (1.0/deltaCoeffs[curGlobalFace]) - ngbDn;

      scalar ownMu = mu.internalField()[curOwner];
      scalar ngbMu = mu.internalField()[curNeighbour];

      scalar ownLambda = lambda.internalField()[curOwner];
      scalar ngbLambda = lambda.internalField()[curNeighbour];

      scalar ownK = (2*ownMu + ownLambda);
      scalar ngbK = (2*ngbMu + ngbLambda);

      // Interface displacement
      // no decomposition of traction - philipc

      // explicit terms are lumped together in Q
      vector Qa =
    ownMu*(ownN&ownGradU.T())
    + (ownLambda*ownTrGradU*ownN)
    - (ownMu+ownLambda)*(ownN&ownGradU);
      vector Qb =
    ngbMu*(ownN&ngbGradU.T())
    + (ngbLambda*ngbTrGradU*ownN)
    - (ngbMu+ngbLambda)*(ownN&ngbGradU);

      // non-orthogonality correction terms

      // delta vectors
      // vector delta = CI[curNeighbour] - CI[curOwner];
      // vector ownDelta = w[curGlobalFace]*delta;
      // vector ngbDelta = delta - ownDelta;
      // it is more accurate to use the face to cell vectors
      // because cells either side of the interface have different
      // non-orthogonality
      vector ownDelta = CfI[curGlobalFace] - CI[curOwner];
      vector ngbDelta = CI[curNeighbour] - CfI[curGlobalFace];

      // correction vectors
      vector ownKcorr = pos(ownN & ownDelta)*(ownN - (ownDelta/ownDn));
      vector ngbKcorr = pos(ownN & ngbDelta)*(ownN - (ngbDelta/ngbDn));
      vector ownCorrection = ownKcorr & ownGradU;
      vector ngbCorrection = ngbKcorr & ngbGradU;
      Qa += (2.0*ownMu + ownLambda) * ( ownCorrection );
      Qb += (2.0*ngbMu + ngbLambda) * ( ngbCorrection );

      vector ownU = UI[curOwner];
      vector ngbU = UI[curNeighbour];

      interU[faceI] =
    (
     ownK*ngbDn*ownU + ngbK*ownDn*ngbU
     + (ownDn*ngbDn*(Qb - Qa))
     )
    /(ownK*ngbDn + ngbK*ownDn);


      // Implicit coupling

      scalar wRevLin = 1.0 - w[curGlobalFace];

      scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);

      scalar Dnf = 1.0/deltaCoeffs[curGlobalFace];

      // Owner
      diag[curOwner] += Kf*magS/Dnf;

      upper[curGlobalFace] -= Kf*magS/Dnf;

      // philipc - no decomposition
      source[curOwner] +=
    (
     ownK*ngbDn*Qb
     + ngbK*ownDn*Qa
     )*magS
    /(ownK*ngbDn + ngbK*ownDn);

      // Neighbour
      diag[curNeighbour] += Kf*magS/Dnf;

      // philipc - no decomposition
      source[curNeighbour] -=
    (
     ownK*ngbDn*Qb
     + ngbK*ownDn*Qa
     )*magS
    /(ownK*ngbDn + ngbK*ownDn);
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
      const vectorField& curProcCfI =
    mesh.Cf().boundaryField()[curPatch];

      const scalarField curProcOwnMu =
    mu.boundaryField()[curPatch].patchInternalField();
      const scalarField curProcNgbMu =
    mu.boundaryField()[curPatch].patchNeighbourField();

      const scalarField curProcOwnLambda =
    lambda.boundaryField()[curPatch].patchInternalField();
      const scalarField curProcNgbLambda =
    lambda.boundaryField()[curPatch].patchNeighbourField();

      const unallocLabelList& curProcFaceCells =
    mesh.boundary()[curPatch].faceCells();

      forAll(processorPatchFaces()[patchI], faceI)
        {
      label curFace = processorPatchFaces()[patchI][faceI];

      scalar ngbDn = curProcW[curFace]*(1.0/curProcDeltaCoeffs[curFace]);
      scalar ownDn = (1.0/curProcDeltaCoeffs[curFace]) - ngbDn;

      scalar magS = curProcMagS[curFace];

      scalar ownMu  = curProcOwnMu[curFace];
      scalar ngbMu  = curProcNgbMu[curFace];

      scalar ownLambda  = curProcOwnLambda[curFace];
      scalar ngbLambda  = curProcNgbLambda[curFace];

      vector ownN = curProcS[curFace]/curProcMagS[curFace];
      //vector ngbN = -ownN;

      // extrapolate tangential gradient
      tensor ownGradU = curProcOwnGradU[curFace];
      tensor ngbGradU = curProcNgbGradU[curFace];
      //tensor ownSGradU = ((I-ownN*ownN)&ownGradU);
      //tensor ngbSGradU = ((I-ngbN*ngbN)&ngbGradU);

      scalar ownTrGradU = tr(ownGradU);
      scalar ngbTrGradU = tr(ngbGradU);


      // Interface displacement

      scalar ownK = (2*ownMu + ownLambda);
      scalar ngbK = (2*ngbMu + ngbLambda);

      // explicit terms
      vector Qa =
        ownMu*(ownN&ownGradU.T())
        + (ownLambda*ownTrGradU*ownN)
        - (ownMu+ownLambda)*(ownN&ownGradU);
      vector Qb =
        ngbMu*(ownN&ngbGradU.T())
        + (ngbLambda*ngbTrGradU*ownN)
        - (ngbMu+ngbLambda)*(ownN&ngbGradU);

      // non-orthogonal correction

      // delta vectors
      //vector delta = curProcNgbC[curFace] - curProcOwnC[curFace];
      //vector ownDelta = curProcW[curFace]*delta;
      //vector ngbDelta = delta - ownDelta;
      // it is more accurate to use the face to cell vectors
      // because cells either side of the interface have different
      // non-orthogonality
      vector ownDelta = curProcCfI[curFace] - curProcOwnC[curFace];
      vector ngbDelta = curProcNgbC[curFace] - curProcCfI[curFace];

      // correction vectors
      vector ownKcorr = pos(ownN & ownDelta)*(ownN - (ownDelta/ownDn));
      vector ngbKcorr = pos(ownN & ngbDelta)*(ownN - (ngbDelta/ngbDn));
      vector ownCorrection = ownKcorr & ownGradU;
      vector ngbCorrection = ngbKcorr & ngbGradU;
      Qa += (2.0*ownMu + ownLambda) * ( ownCorrection );
      Qb += (2.0*ngbMu + ngbLambda) * ( ngbCorrection );

      vector ownU = curProcOwnU[curFace];
      vector ngbU = curProcNgbU[curFace];

      curProcInterU[faceI] =
        (
         ownK*ngbDn*ownU + ngbK*ownDn*ngbU
         + (ownDn*ngbDn*(Qb - Qa))
         )
        /(ownK*ngbDn + ngbK*ownDn);


      // Implicit coupling

      scalar wRevLin = 1.0 - curProcW[curFace];

      scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);

      scalar Dnf = 1.0/curProcDeltaCoeffs[curFace];

      // Owner
      diag[curProcFaceCells[curFace]] += Kf*magS/Dnf;

      boundaryCoeffs[curPatch][curFace] += Kf*magS*vector::one/Dnf;

      // no decomposition
      source[curProcFaceCells[curFace]] +=
        (
         ownK*ngbDn*Qb
         + ngbK*ownDn*Qa
         )*magS
        /(ownK*ngbDn + ngbK*ownDn);
        }
    }
}

// ************************************************************************* //
