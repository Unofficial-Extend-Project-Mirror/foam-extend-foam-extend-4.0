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
    Calculates the force and average displacement of the specified patch
    and writes it to a file

Author
    philip.cardiff@ucd.ie

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  Foam::argList::validOptions.insert("noMeshUpdate", "");
  Foam::argList::validOptions.insert("nonLinear", "");
  argList::validArgs.append("patchName");

# include "addTimeOptions.H"

# include "setRootCase.H"

# include "createTime.H"

  word patchName(args.additionalArgs()[0]);

  // Get times list
  instantList Times = runTime.times();

  // set startTime and endTime depending on -time and -latestTime options
# include "checkTimeOptions.H"

  runTime.setTime(Times[startTime], startTime);

# include "createMesh.H"

  bool noMeshUpdate = args.optionFound("noMeshUpdate");
  bool nonLinear = args.optionFound("nonLinear");

  // check patch exists
  label patchID = mesh.boundaryMesh().findPatchID(patchName);
  if (patchID == -1)
    {
      FatalError
          << "Cannot find patch " << patchName << exit(FatalError);
    }

  // open file
  OFstream forceDispFile("forceDisp_"+patchName+".dat");
  label width = 20;
  forceDispFile << "average disp";
  forceDispFile.width(width);
  forceDispFile << "totalForce";
  forceDispFile << endl;

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

      IOobject Uheader
          (
              "U",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ
              );

      // Check sigma exists
      if (sigmaheader.headerOk() && Uheader.headerOk())
      {
          Info<< "    Reading sigma" << endl;
          volSymmTensorField sigma(sigmaheader, mesh);
          Info<< "    Reading U" << endl;
          volVectorField U(Uheader, mesh);

          Info << nl;

          // calculate patch force
          const vectorField n = mesh.boundary()[patchID].nf();
          const vectorField& Sf = mesh.boundary()[patchID].Sf();
          const symmTensorField& sigmaPatch = sigma.boundaryField()[patchID];
          vector totalForce = vector::zero;
          scalar normalForce = 0.0;
          scalar shearForce = 0.0;
          if (nonLinear)
          {
              // assuming sigma is the second Piola-Kirchhoff tensor

              // get deformation gradient
              volTensorField gradU
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
              tensorField F = I + gradU.boundaryField()[patchID];
              vectorField force = Sf & (sigmaPatch & F);
              tensorField Finv = inv(F);
              vectorField nCur = Finv & n;
              nCur /= mag(nCur);
              normalForce = sum(nCur & force);
              shearForce = sum( mag((I - sqr(nCur)) & force) );
              totalForce = sum(force);
          }
          else
          {
              // small strain or UL large strain
              // (as the mesh is moved and sigma is Cauchy)
              //Info << "\tSmall Strain Total Lagrangian"<<nl<<endl;
              vectorField force = Sf & sigmaPatch;
              normalForce = sum(n & force);
              shearForce = sum( mag((I - sqr(n)) & force) );
              totalForce = sum(force);
          }

          // calculate average displacement
          // these should be a weighted average but OK for most models
          vector disp = average(U.boundaryField()[patchID]);

          // write to file
          forceDispFile
              << disp.x() << " " << disp.y() << " " << disp.z();
          forceDispFile.width(width);
          forceDispFile
              << totalForce.x() << " " << totalForce.y()
              << " " << totalForce.z();
          forceDispFile
              << " " << normalForce << " " << shearForce
              << endl;
      }
  }

  Info << nl << "End" << endl;

  return 0;
}


// ************************************************************************* //
