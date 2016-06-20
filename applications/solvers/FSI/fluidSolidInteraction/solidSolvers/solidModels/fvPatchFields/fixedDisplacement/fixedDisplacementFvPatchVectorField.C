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

#include "fixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    secondOrder_(false)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    secondOrder_(ptf.secondOrder_)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    secondOrder_(false)
{
    if (dict.found("secondOrder"))
    {
        secondOrder_ = Switch(dict.lookup("secondOrder"));
        Info << "Second order correction: " << secondOrder_ << endl;
    }

    Info << "Creating fixed displacement boundary condition" << endl;
}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    secondOrder_(pivpvf.secondOrder_)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    secondOrder_(pivpvf.secondOrder_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> > fixedDisplacementFvPatchVectorField::
snGrad() const
{
    word DName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    if (secondOrder_)
    {
        vectorField dDP = (k&gradD.patchInternalField());
        vectorField nGradDP = (n&gradD.patchInternalField());

        return
            2
           *(
                *this
              - (this->patchInternalField() + dDP)
            )*this->patch().deltaCoeffs()
          - nGradDP;
    }

    return
    (
        *this
      - (patchInternalField() + (k&gradD.patchInternalField()))
    )*this->patch().deltaCoeffs();
}

tmp<Field<vector> > fixedDisplacementFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    word DName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    if (secondOrder_)
    {
        vectorField dDP = (k&gradD.patchInternalField());
        vectorField nGradDP = (n&gradD.patchInternalField());

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dDP)
              - this->patchInternalField()
            )
          - nGradDP;
    }

    return this->patch().deltaCoeffs()
       *(*this - (k&gradD.patchInternalField()));
}

void fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
