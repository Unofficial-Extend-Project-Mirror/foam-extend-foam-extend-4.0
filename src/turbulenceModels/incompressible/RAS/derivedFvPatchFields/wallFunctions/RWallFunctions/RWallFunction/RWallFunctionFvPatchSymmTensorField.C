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

\*---------------------------------------------------------------------------*/

#include "RWallFunctionFvPatchSymmTensorField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void RWallFunctionFvPatchSymmTensorField::checkType()
{
    if (!this->patch().isWall())
    {
        FatalErrorIn("RWallFunctionFvPatchSymmTensorField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void RWallFunctionFvPatchSymmTensorField::setInInternalField
(
    const symmTensorField& f
) const
{
    symmTensorField& vf = const_cast<symmTensorField& >
    (
        this->internalField()
    );

    const labelList& faceCells = this->patch().faceCells();

    // Apply the refValue into the cells next to the boundary
    forAll (faceCells, faceI)
    {
        vf[faceCells[faceI]] = f[faceI];
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RWallFunctionFvPatchSymmTensorField::RWallFunctionFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    zeroGradientFvPatchField<symmTensor>(p, iF),
    UName_("U"),
    nutName_("nut")
{
    checkType();
}


RWallFunctionFvPatchSymmTensorField::RWallFunctionFvPatchSymmTensorField
(
    const RWallFunctionFvPatchSymmTensorField& ptf,
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<symmTensor>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    nutName_(ptf.nutName_)
{
    checkType();
}


RWallFunctionFvPatchSymmTensorField::RWallFunctionFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<symmTensor>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut"))
{
    checkType();
}


RWallFunctionFvPatchSymmTensorField::RWallFunctionFvPatchSymmTensorField
(
    const RWallFunctionFvPatchSymmTensorField& Rwfpsf
)
:
    zeroGradientFvPatchField<symmTensor>(Rwfpsf),
    UName_(Rwfpsf.UName_),
    nutName_(Rwfpsf.nutName_)
{
    checkType();
}


RWallFunctionFvPatchSymmTensorField::RWallFunctionFvPatchSymmTensorField
(
    const RWallFunctionFvPatchSymmTensorField& Rwfpsf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    zeroGradientFvPatchField<symmTensor>(Rwfpsf, iF),
    UName_(Rwfpsf.UName_),
    nutName_(Rwfpsf.nutName_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void RWallFunctionFvPatchSymmTensorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    const scalarField& nutw =
        lookupPatchField<volScalarField, scalar>(nutName_);

    const fvPatchVectorField& Uw =
        lookupPatchField<volVectorField, vector>(UName_);

    // Old treatment: anisotropy not enforced
    vectorField snGradU = Uw.snGrad();

    tensorField gradUw = patch().nf()*snGradU;

    // Calculate near-wall shear-stress tensor
    symmTensorField tauw = -nutw*2*symm(gradUw);

//     // Normal coordinate
//     vectorField n = patch().nf();

//     // Tangential coordinate.  Take care of mag(U) = 0;
//     vectorField Utan = (I - sqr(n)) & Uw.patchInternalField();
//     scalarField magUtan = mag(Utan);
//     Utan /= magUtan + SMALL;

//     vectorField t = pos(magUtan - SMALL)*Utan
//         + neg(magUtan - SMALL)*vector(1, 0, 0);

//     // Binormal coordinate
//     vectorField l = n ^ t;

    // Reset the shear components of the stress tensor
    replace(symmTensor::XY, tauw.component(symmTensor::XY));
    replace(symmTensor::XZ, tauw.component(symmTensor::XZ));
    replace(symmTensor::YZ, tauw.component(symmTensor::YZ));
}


void RWallFunctionFvPatchSymmTensorField::write(Ostream& os) const
{
    zeroGradientFvPatchSymmTensorField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchSymmTensorField,
    RWallFunctionFvPatchSymmTensorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
