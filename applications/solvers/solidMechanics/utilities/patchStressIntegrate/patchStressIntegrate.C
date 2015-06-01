/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    Calculates the total forces on a patch:
            total force vector
            total normal force
            total force in each direction (x, y and z)

Author
    philip.cardiff@ucd.ie

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  argList::validArgs.append("patch name");
  argList::validOptions.insert("noMeshUpdate", "");
  argList::validOptions.insert("nonLinear", "");

# include "addTimeOptions.H"
# include "setRootCase.H"
# include "createTime.H"

  // Get times list
  instantList Times = runTime.times();

  // set startTime and endTime depending on -time and -latestTime options
# include "checkTimeOptions.H"

  runTime.setTime(Times[startTime], startTime);

# include "createMesh.H"

  bool noMeshUpdate = args.optionFound("noMeshUpdate");
  bool nonLinear = args.optionFound("nonLinear");

  word patchName(args.additionalArgs()[0]);
  label patchID = mesh.boundaryMesh().findPatchID(patchName);
  if (patchID == -1)
  {
      FatalError << "Cannot find patch " << patchName
                 << exit(FatalError);
  }

  for (label i=startTime; i<endTime; i++)
  {
      runTime.setTime(Times[i], i);

      Info<< "Time = " << runTime.timeName() << endl;

      if (!noMeshUpdate)
      {
          mesh.readUpdate();
      }

      IOobject sigmaheader
          (
              "sigma",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ
              );

      // Check sigma exists
      if (sigmaheader.headerOk())
      {
          Info<< "\tReading sigma" << endl;
          volSymmTensorField sigma(sigmaheader, mesh);

          Info<< nl;

          // gradU needed for nonLinear
          volTensorField* gradUPtr = NULL;
          volSymmTensorField* sigmaCauchyPtr = NULL;
          if (nonLinear)
          {
              gradUPtr = new volTensorField
                  (
                      IOobject
                      (
                          "grad(U)",
                          runTime.timeName(),
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::NO_WRITE
                          ),
                      mesh
                      );

              sigmaCauchyPtr = new volSymmTensorField
                  (
                      IOobject
                      (
                          "sigmaCauchy",
                          runTime.timeName(),
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::NO_WRITE
                          ),
                      mesh
                      );
          }

          //vector netForce = vector::zero;
          //vector netForceCauchy = vector::zero;
          //scalar maxPatchForce = 0.0;
          //forAll(mesh.boundary(), patchID)
          {
              vectorField n = mesh.boundary()[patchID].nf();
              const vectorField& Sf = mesh.boundary()[patchID].Sf();
              const symmTensorField& sigmaPatch =
                  sigma.boundaryField()[patchID];

              vectorField totalForce(sigmaPatch.size(), vector::zero);
              scalar totalNormalForce = 0.0;
              vector totalShearForce = vector::zero;

              vectorField totalForceCauchy(sigmaPatch.size(), vector::zero);
              scalar totalNormalForceCauchy = 0.0;
              vector totalShearForceCauchy = vector::zero;
              if (nonLinear)
              {
                  // Note: only for TL models, not correct for UL
                  // models yet - todo

                  // We use two separate methods to calculate the force
                  // for the nonlinear models
                  // both methods should be equivalent
                  // they are both used just to check everything is as it
                  // should be
                  // Force == currentAreas & sigmaCauchy ==
                  // referenceArea & sigma2PK & deformationGradient

                  // deformation gradient
                  tensorField F = I + gradUPtr->boundaryField()[patchID];

                  const scalarField J = det(F);
                  const tensorField Finv = hinv(F);
                  // current deformed patch area vectors are given by
                  // Nanson's formula
                  const vectorField deformedSf = J * Finv & Sf;
                  const vectorField deformedN = deformedSf/mag(deformedSf);
                  const symmTensorField& sigmaCauchyPatch =
                      sigmaCauchyPtr->boundaryField()[patchID];

                  // reference areas and 2nd Piola-Kirchhoff stress
                  totalForce = Sf & (sigmaPatch & F);
                  totalNormalForce = sum(deformedN & (totalForce));
                  totalShearForce = sum((I -sqr(deformedN)) & (totalForce));

                  // deformed normals and Cauchy stress
                  totalForceCauchy = deformedSf & (sigmaCauchyPatch);
                  totalNormalForceCauchy = sum(deformedN & (totalForceCauchy));
                  totalShearForceCauchy =
                      sum((I -sqr(deformedN)) & (totalForceCauchy));
                  //netForceCauchy += sum(totalForceCauchy);
              }
              else
              {
                  // small strain
                  totalForce = Sf & sigmaPatch;
                  totalNormalForce = sum(n & (totalForce));
                  totalShearForce = sum((I -sqr(n)) & (totalForce));
              }
              //netForce += sum(totalForce);
              // scalar totalNormalForce = sum(n & (totalForce));
              // vector totalShearForce = sum((I -sqr(n)) & (totalForce));

              //maxPatchForce = max(maxPatchForce, mag(sum(totalForce)));

              Info<< "Patch: " << mesh.boundary()[patchID].name() << nl
                  << "\tTotal Force:\t\t" << sum(totalForce) << " N\n"
                  << "\tTotal Normal Force:\t" << totalNormalForce <<  " N\n"
                  << "\tTotal Shear Force:\t" << totalShearForce <<  " N\n";
              if (nonLinear)
              {
                  Info<< "\tForces calculated with Cauchy stress\n"
                      << "\tTotal Force:\t\t" << sum(totalForceCauchy) << " N\n"
                      << "\tTotal Normal Force:\t" << totalNormalForceCauchy
                      <<  " N\n"
                      << "\tTotal Shear Force:\t" << totalShearForceCauchy
                      <<  " N\n";
              }
              Info<< endl;
          }
          // scalar percentNetForce = 100.0*mag(netForce)/maxPatchForce;
          // scalar percentNetForceCauchy =
          // 100.0*mag(netForceCauchy)/maxPatchForce;
          // Info<< nl << "Net force on model is "
          // << netForce << " N\twhich is "
          //       << percentNetForce << "% of maximum patch force";

          Info<< endl;
      }
  }

  Info<< nl << "End" << endl;

  return 0;
}


// ************************************************************************* //
