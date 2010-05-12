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

#include "BICCG.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BICCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<BICCG>
        addBICCGSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<BICCG>
        addBICCGAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IStringStream Foam::BICCG::solverDataStream
(
    const scalar tolerance,
    const scalar relTol
) const
{
    return IStringStream
    (
        "{ preconditioner { type DILU; }"
        "  minIter 0; maxIter 1000; "
        "  tolerance " + name(tolerance) + "; relTol " + name(relTol) + "; }"
    );
}


Foam::IStringStream Foam::BICCG::solverDataStream
(
    const word& solverName,
    Istream& solverData
) const
{
    scalar tolerance(readScalar(solverData));
    scalar relTol(readScalar(solverData));

    return IStringStream
    (
        solverName + "{ preconditioner { type DILU; }"
        "  minIter 0; maxIter 1000; "
        "  tolerance " + name(tolerance) + "; relTol " + name(relTol) + "; }"
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BICCG::BICCG
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
    PBiCG
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        solverDataStream(tolerance, relTol)()
    )
{}


Foam::BICCG::BICCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    Istream& solverData
)
:
    PBiCG
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

void Foam::BICCG::read(Istream& solverData)
{
    word solverName(solverData);
    PBiCG::read(solverDataStream(solverName, solverData)());
}


// ************************************************************************* //
