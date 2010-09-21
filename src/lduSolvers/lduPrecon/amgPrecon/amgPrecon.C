/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    amgPrecon

Description
    AMG preconditioner

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "amgPrecon.H"
#include "fineAmgLevel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(amgPrecon, 0);

    lduPreconditioner::
        addsymMatrixConstructorToTable<amgPrecon>
        addamgPreconditionerSymMatrixConstructorToTable_;

    lduPreconditioner::
        addasymMatrixConstructorToTable<amgPrecon>
        addamgPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::amgPrecon::amgPrecon
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaceFields,
    const dictionary& dict
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaceFields
    ),
    cycle_(amgCycle::cycleNames_.read(dict.lookup("cycle"))),
    nPreSweeps_(readLabel(dict.lookup("nPreSweeps"))),
    nPostSweeps_(readLabel(dict.lookup("nPostSweeps"))),
    nMaxLevels_(readLabel(dict.lookup("nMaxLevels"))),
    scale_(dict.lookup("scale")),
    amgPtr_
    (
        new amgCycle
        (
            autoPtr<amgLevel>
            (
                new fineAmgLevel
                (
                    matrix,
                    coupleBouCoeffs,
                    coupleIntCoeffs,
                    interfaceFields,
                    dict,
                    dict.lookup("policy"),
                    readLabel(dict.lookup("groupSize")),
                    readLabel(dict.lookup("minCoarseEqns")),
                    dict.lookup("smoother")
                )
            )
        )
    ),
    xBuffer_(matrix.lduAddr().size())
{
    // Make coarse levels
    amgPtr_->makeCoarseLevels(nMaxLevels_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::amgPrecon::~amgPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::amgPrecon::nLevels() const
{
    return amgPtr_->nLevels();
}


const Foam::scalarField& Foam::amgPrecon::residual
(
    const scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Calculate residual
    amgPtr_->residual(x, b, cmpt, xBuffer_);

    return xBuffer_;
}


void Foam::amgPrecon::cycle
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    amgPtr_->fixedCycle
    (
        x,
        b,
        cmpt,
        xBuffer_,
        cycle_,
        nPreSweeps_,
        nPostSweeps_,
        scale_
    );
}


void Foam::amgPrecon::precondition
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Execute preconditioning
    residual(x, b, cmpt);
    cycle(x, b, cmpt);
}


// ************************************************************************* //
