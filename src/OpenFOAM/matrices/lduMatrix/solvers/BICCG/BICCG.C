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

    lduSolver::addsymMatrixConstructorToTable<BICCG>
        addBICCGSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<BICCG>
        addBICCGAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::dictionary Foam::BICCG::solverDict
(
    const scalar tol,
    const scalar relTol
)
{
    dictionary dict(IStringStream("solver PBiCG; preconditioner DILU;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


Foam::dictionary Foam::BICCG::solverDict
(
    Istream& is
)
{
    scalar tol(readScalar(is));
    scalar relTol(readScalar(is));

    return solverDict(tol, relTol);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BICCG::BICCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    PBiCG
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    )
{}


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
        solverDict(tolerance, relTol)
    )
{}

// ************************************************************************* //
