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
    liquidFilmFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids in
    liquid film formulation.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFaMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "readTransportProperties.H"
#   include "createFaFields.H"
#   include "createFvFields.H"

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readSolutionControls.H"
#       include "readTimeControls.H"
#       include "surfaceCourantNo.H"
#       include "capillaryCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int iCorr = 0; iCorr < nCorr; iCorr++)
        {
            phi2s = fac::interpolate(h)*phis;

#           include "calcFrictionFactor.H"

            faVectorMatrix UsEqn
            (
                fam::ddt(h, Us)
              + fam::div(phi2s, Us)
              + fam::Sp(0.0125*frictionFactor*mag(Us), Us)
             ==
                Gs*h
              - fam::Sp(Sd, Us)
            );

            UsEqn.relax();
            solve(UsEqn == - fac::grad(ps*h)/rhol + ps*fac::grad(h)/rhol);

            areaScalarField UsA = UsEqn.A();

            Us = UsEqn.H()/UsA;
            Us.correctBoundaryConditions();

            phis = (fac::interpolate(Us) & aMesh.Le())
              - fac::interpolate(1.0/(rhol*UsA))
               *fac::lnGrad(ps*h)*aMesh.magLe()
              + fac::interpolate(ps/(rhol*UsA))
               *fac::lnGrad(h)*aMesh.magLe();

            faScalarMatrix hEqn
            (
                fam::ddt(h)
              + fam::div(phis, h)
             ==
                Sm
              - fam::Sp
                (
                    Sd/(h + dimensionedScalar("small", dimLength, SMALL)),
                    h
                )
            );

            hEqn.relax();
            hEqn.solve();

            phi2s = hEqn.flux();

            // Bound h
            h.internalField() = max
            (
                max
                (
                    h.internalField(),
                    fac::average(max(h, h0))().internalField()
                    *pos(h0.value() - h.internalField())
                ),
                h0.value()
            );

            ps = rhol*Gn*h - sigma*fac::laplacian(h);
            ps.correctBoundaryConditions();

            Us -= (1.0/(rhol*UsA))*fac::grad(ps*h)
                - (ps/(rhol*UsA))*fac::grad(h);
            Us.correctBoundaryConditions();
        }

        if (runTime.outputTime())
        {
            vsm.mapToVolume(h, H.boundaryField());
            vsm.mapToVolume(Us, U.boundaryField());

            runTime.write();
        }

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    return(0);
}

// ************************************************************************* //
