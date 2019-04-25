/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar velMag = 0.0;

// HR 26.06.18: A parallel run has at least two cells and therefore at least
// one internal face in the global mesh. It may be a processor boundary, but
// this is captured by max(mag(phi)).
// Old formulation hangs on parallel cases where one partition is degenerated
// to a single cell.
if (mesh.nInternalFaces() || Pstream::parRun())
{
    surfaceScalarField magPhi = mag(phi);

    surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*magPhi;

    const scalar deltaT = runTime.deltaT().value();

    CoNum = max(SfUfbyDelta/mesh.magSf()).value()*deltaT;

    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf())).value()*deltaT;

    velMag = max(magPhi/mesh.magSf()).value();
}
else
{
    // Single cell mesh: Co is still defined; use cell formulation

    const scalar deltaT = runTime.deltaT().value();

    const scalar deltaX = Foam::cbrt(mesh.V()[0]);

    // recover velocity field in a more general way
    const volVectorField& URef
	    = mesh.db().lookupObject<const volVectorField>("U");

    velMag = mag(URef[0]);

    CoNum = velMag*deltaT/deltaX;

    meanCoNum = CoNum;
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum
    << " velocity magnitude: " << velMag
    << endl;

// ************************************************************************* //
