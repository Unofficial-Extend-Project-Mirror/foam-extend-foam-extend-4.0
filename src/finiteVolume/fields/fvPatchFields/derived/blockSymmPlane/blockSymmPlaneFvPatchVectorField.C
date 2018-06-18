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

\*---------------------------------------------------------------------------*/

#include "blockSymmPlaneFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockSymmPlaneFvPatchVectorField::blockSymmPlaneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    basicSymmetryFvPatchVectorField(p, iF)
{}


Foam::blockSymmPlaneFvPatchVectorField::blockSymmPlaneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchVectorField(p, iF, dict)
{}


Foam::blockSymmPlaneFvPatchVectorField::blockSymmPlaneFvPatchVectorField
(
    const blockSymmPlaneFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::blockSymmPlaneFvPatchVectorField::blockSymmPlaneFvPatchVectorField
(
    const blockSymmPlaneFvPatchVectorField& ptf
)
:
    basicSymmetryFvPatchVectorField(ptf)
{}


Foam::blockSymmPlaneFvPatchVectorField::blockSymmPlaneFvPatchVectorField
(
    const blockSymmPlaneFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    basicSymmetryFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Disabled: see coupled coefficient
Foam::tmp<Foam::vectorField>
Foam::blockSymmPlaneFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<vectorField>(new vectorField(this->size(), vector::zero));
}


// Disabled: see coupled coefficient
Foam::tmp<Foam::vectorField>
Foam::blockSymmPlaneFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<vectorField>(new vectorField(this->size(), vector::zero));
}


// Disabled: see coupled coefficient
Foam::tmp<Foam::vectorField>
Foam::blockSymmPlaneFvPatchVectorField::gradientInternalCoeffs() const
{
    return tmp<vectorField>(new vectorField(this->size(), vector::zero));
}


// Disabled: see coupled coefficient
Foam::tmp<Foam::vectorField>
Foam::blockSymmPlaneFvPatchVectorField::gradientBoundaryCoeffs() const
{
    return tmp<vectorField>(new vectorField(this->size(), vector::zero));
}


// Coupled coefficient
Foam::tmp<Foam::vectorCoeffField>
Foam::blockSymmPlaneFvPatchVectorField::blockValueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<vectorCoeffField> tcoeff
    (
        new vectorCoeffField(this->size())
    );

    tcoeff().asSquare() = I - this->patch().nf()*this->patch().nf();

    return tcoeff;
}


// Coupled coefficient
Foam::tmp<Foam::vectorField>
Foam::blockSymmPlaneFvPatchVectorField::blockValueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        *this
      - transform
        (
            I - sqr(this->patch().nf()),
            this->patchInternalField()
        );
}


// Coupled coefficient
Foam::tmp<Foam::vectorCoeffField>
Foam::blockSymmPlaneFvPatchVectorField::blockGradientInternalCoeffs() const
{
    tmp<vectorCoeffField> tcoeff
    (
        new vectorCoeffField(this->size())
    );

    tcoeff().asSquare() =
        -this->patch().deltaCoeffs()*this->patch().nf()*this->patch().nf();

    return tcoeff;
}


// Coupled coefficient
Foam::tmp<Foam::vectorField>
Foam::blockSymmPlaneFvPatchVectorField::blockGradientBoundaryCoeffs() const
{
    return
    (
        sqr(this->patch().nf())
      & (
            this->snGrad()
          - transform
            (
                -this->patch().deltaCoeffs()*sqr(this->patch().nf()),
                this->patchInternalField()
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    blockSymmPlaneFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
