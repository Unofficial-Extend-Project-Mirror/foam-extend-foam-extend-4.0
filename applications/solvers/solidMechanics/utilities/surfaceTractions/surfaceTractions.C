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

Application
    surfaceTractions

Description
    Calculates and writes the surface tractions as a volVectorField, using
    the sigma volSymmTensorField

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  Foam::argList::validOptions.insert("nonLinear", "");

# include "addTimeOptions.H"
# include "setRootCase.H"
# include "createTime.H"

  bool nonLinear = args.optionFound("nonLinear");

  // Get times list
  instantList Times = runTime.times();

  // set startTime and endTime depending on -time and -latestTime options
# include "checkTimeOptions.H"

  runTime.setTime(Times[startTime], startTime);

# include "createMesh.H"

  for (label i=startTime; i<endTime; i++)
  {
      runTime.setTime(Times[i], i);

      Info<< "Time = " << runTime.timeName() << endl;

      mesh.readUpdate();

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
          mesh.readUpdate();

          Info<< "    Reading sigma" << endl;
          volSymmTensorField sigma(sigmaheader, mesh);

          surfaceVectorField n = mesh.Sf()/mesh.magSf();

          volVectorField totalTraction
              (
                  IOobject
                  (
                      "totalTraction",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                      ),
                  mesh,
                  dimensionedVector("zero", dimForce/dimArea, vector::zero)
                  );
          volScalarField normalTraction
              (
                  IOobject
                  (
                      "normalTraction",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                      ),
                  mesh,
                  dimensionedScalar("zero", dimForce/dimArea, 0.0)
                  );
          volVectorField shearTraction
              (
                  IOobject
                  (
                      "shearTraction",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                      ),
                  mesh,
                  dimensionedVector("zero", dimForce/dimArea, vector::zero)
                  );

          volTensorField* gradUPtr = NULL;
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
          }

          forAll(totalTraction.boundaryField(), patchi)
          {
              const vectorField& nb = n.boundaryField()[patchi];
              const symmTensorField& sigmab = sigma.boundaryField()[patchi];

              if (nonLinear)
              {
                  tensorField F = I + gradUPtr->boundaryField()[patchi];
                  totalTraction.boundaryField()[patchi] = nb & (sigmab & F);
              }
              else
              {
                  totalTraction.boundaryField()[patchi] = nb & sigmab;
              }
              normalTraction.boundaryField()[patchi] =
                  nb & totalTraction.boundaryField()[patchi];
              shearTraction.boundaryField()[patchi] =
                  (I -sqr(nb)) & totalTraction.boundaryField()[patchi];
          }
          totalTraction.write();
          normalTraction.write();
          shearTraction.write();
      }
      else
      {
          Info<< "    No sigma field" << endl;
      }
      Info<< endl;
  }

  Info<< "End" << endl;

  return(0);
}


// ************************************************************************* //
