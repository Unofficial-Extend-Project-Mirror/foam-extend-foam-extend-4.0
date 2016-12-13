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

#include "GAMGPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGPreconditioner, 0);

    lduPreconditioner::addsymMatrixConstructorToTable
    <GAMGPreconditioner> addGAMGPreconditionerSymMatrixConstructorToTable_;

    lduPreconditioner::addasymMatrixConstructorToTable
    <GAMGPreconditioner> addGAMGPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGPreconditioner::GAMGPreconditioner
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    GAMG_
    (
        "dummy",
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    ),
    nVcycles_(2)
{
    readControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGPreconditioner::~GAMGPreconditioner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GAMGPreconditioner::readControls()
{
    GAMG_.readControls();
    nVcycles_ = GAMG_.dict().lookupOrDefault<label>("nVcycles", 2);
}


void Foam::GAMGPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction cmpt
) const
{
    wA = 0.0;
    scalarField AwA(wA.size());
    scalarField finestCorrection(wA.size());
    scalarField finestResidual(rA);

    // Create coarse grid correction fields
    PtrList<scalarField> coarseCorrX;

    // Create coarse grid sources
    PtrList<scalarField> coarseB;

    // Create the smoothers for all levels
    PtrList<lduSmoother> smoothers;

    // Initialise the above data structures
    GAMG_.initVcycle(coarseCorrX, coarseB, smoothers);

    for (label cycle=0; cycle<nVcycles_; cycle++)
    {
        GAMG_.Vcycle
        (
            smoothers,
            wA,
            rA,
            AwA,
            finestCorrection,
            finestResidual,
            coarseCorrX,
            coarseB,
            cmpt
        );

        if (cycle < nVcycles_ - 1)
        {
            // Calculate finest level residual field
            matrix_.Amul(AwA, wA, coupleBouCoeffs_, interfaces_, cmpt);
            finestResidual = rA;
            finestResidual -= AwA;
        }
    }
}


// ************************************************************************* //
