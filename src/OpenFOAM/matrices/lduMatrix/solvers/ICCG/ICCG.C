/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "ICCG.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ICCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<ICCG>
        addICCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IStringStream Foam::ICCG::solverDataStream
(
    const scalar tolerance,
    const scalar relTol
) const
{
    return IStringStream
    (
        "{ preconditioner { type DIC; }"
        "  minIter 0; maxIter 5000; "
        "  tolerance " + name(tolerance) + "; relTol " + name(relTol) + "; }"
    );
}


Foam::IStringStream Foam::ICCG::solverDataStream
(
    const word& solverName,
    Istream& solverData
) const
{
    scalar tolerance(readScalar(solverData));
    scalar relTol(readScalar(solverData));

    return IStringStream
    (
        solverName + "{ preconditioner { type DIC; }"
        "  minIter 0; maxIter 5000; "
        "  tolerance " + name(tolerance) + "; relTol " + name(relTol) + "; }"
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ICCG::ICCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalar tolerance,
    const scalar relTol
)
:
    PCG
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        solverDataStream(tolerance, relTol)()
    )
{}


Foam::ICCG::ICCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    Istream& solverData
)
:
    PCG
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        solverDataStream(word::null, solverData)()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ICCG::read(Istream& solverData)
{
    word solverName(solverData);
    PCG::read(solverDataStream(solverName, solverData)());
}


// ************************************************************************* //
