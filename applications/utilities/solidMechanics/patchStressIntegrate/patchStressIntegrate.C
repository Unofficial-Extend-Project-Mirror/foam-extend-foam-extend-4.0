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
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMesh.H"

   while(runTime.loop())
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;

  volSymmTensorField sigma
    (
     IOobject
     (
      "sigma",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      ),
     mesh
     );

  Info << nl;

  vector netForce = vector::zero;

  forAll(mesh.boundary(), patchID)
    {
      vectorField n = mesh.boundary()[patchID].nf();
      const vectorField& Sf = mesh.boundary()[patchID].Sf();
      const symmTensorField& sigmaPatch = sigma.boundaryField()[patchID];

      vector totalForce = sum(Sf & sigmaPatch);
      netForce += totalForce;
      scalar totalNormalForce = sum(n & (Sf & sigmaPatch));
      vector totalShearForce = sum((I -sqr(n)) & (Sf & sigmaPatch));

      Info << "Patch: " << mesh.boundary()[patchID].name() << nl
          << "\tTotal Force:\t\t" << totalForce << " N\n"
          << "\tTotal Normal Force:\t" << totalNormalForce <<  " N\n"
          << "\tTotal Shear Force:\t" << totalShearForce <<  " N\n" << endl;
    }

  Info << nl << "Net force on model is " << netForce << " N" << endl;

    }

  Info << nl << "End" << endl;

  return 0;
}


// ************************************************************************* //
