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
    blockCoupledScalarTransportFoam

Description
    Solves two coupled transport equations in a block-coupled manner

        1) transport equation for a passive scalar
        2) diffusion only

    This resembles heat exchanging flow through a porous medium

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fieldTypes.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "fvBlockMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix TEqn
            (
                fvm::div(phi, T)
              - fvm::laplacian(DT, T)
             ==
                 alpha*Ts
              - fvm::Sp(alpha, T)
            );

            TEqn.relax();

            fvScalarMatrix TsEqn
            (
              - fvm::laplacian(DTs, Ts)
             ==
                 alpha*T
              - fvm::Sp(alpha, Ts)
            );

            TsEqn.relax();

            // Prepare block system
            fvBlockMatrix<vector2> blockM(blockT);

            // Insert equations into block Matrix
            blockM.insertEquation(0, TEqn);
            blockM.insertEquation(1, TsEqn);

            // Add off-diagonal coupling terms
            scalarField coupling(mesh.nCells(), -alpha.value());

            blockM.insertEquationCoupling(0, 1, coupling);
            blockM.insertEquationCoupling(1, 0, coupling);

            // Update source coupling: coupling terms eliminated from source
            blockM.updateSourceCoupling();

            //- Block coupled solver call
            blockM.solve();

            // Retrieve solution
            blockM.retrieveSolution(0, T.internalField());
            blockM.retrieveSolution(1, Ts.internalField());

            T.correctBoundaryConditions();
            Ts.correctBoundaryConditions();
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
