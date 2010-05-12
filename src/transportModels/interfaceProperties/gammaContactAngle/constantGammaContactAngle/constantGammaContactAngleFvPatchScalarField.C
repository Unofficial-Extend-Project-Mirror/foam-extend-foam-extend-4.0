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

#include "constantGammaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantGammaContactAngleFvPatchScalarField::
constantGammaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    gammaContactAngleFvPatchScalarField(p, iF),
    theta0_(0.0)
{}


constantGammaContactAngleFvPatchScalarField::
constantGammaContactAngleFvPatchScalarField
(
    const constantGammaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    gammaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_)
{}


constantGammaContactAngleFvPatchScalarField::
constantGammaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    gammaContactAngleFvPatchScalarField(p, iF),
    theta0_(readScalar(dict.lookup("theta0")))
{
    evaluate();
}


constantGammaContactAngleFvPatchScalarField::
constantGammaContactAngleFvPatchScalarField
(
    const constantGammaContactAngleFvPatchScalarField& gcpsf
)
:
    gammaContactAngleFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_)
{}


constantGammaContactAngleFvPatchScalarField::
constantGammaContactAngleFvPatchScalarField
(
    const constantGammaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    gammaContactAngleFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> constantGammaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    return tmp<scalarField>(new scalarField(size(), theta0_));
}


void constantGammaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    constantGammaContactAngleFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
