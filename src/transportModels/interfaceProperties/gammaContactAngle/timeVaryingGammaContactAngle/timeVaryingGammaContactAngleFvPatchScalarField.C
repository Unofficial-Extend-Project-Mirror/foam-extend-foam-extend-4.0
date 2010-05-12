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

#include "timeVaryingGammaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingGammaContactAngleFvPatchScalarField::
timeVaryingGammaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    gammaContactAngleFvPatchScalarField(p, iF),
    t0_(0.0),
    thetaT0_(0.0),
    te_(0.0),
    thetaTe_(0.0)
{}


timeVaryingGammaContactAngleFvPatchScalarField::
timeVaryingGammaContactAngleFvPatchScalarField
(
    const timeVaryingGammaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    gammaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    t0_(gcpsf.t0_),
    thetaT0_(gcpsf.thetaT0_),
    te_(gcpsf.te_),
    thetaTe_(gcpsf.te_)
{}


timeVaryingGammaContactAngleFvPatchScalarField::
timeVaryingGammaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    gammaContactAngleFvPatchScalarField(p, iF),
    t0_(readScalar(dict.lookup("t0"))),
    thetaT0_(readScalar(dict.lookup("thetaT0"))),
    te_(readScalar(dict.lookup("te"))),
    thetaTe_(readScalar(dict.lookup("thetaTe")))
{
    evaluate();
}


timeVaryingGammaContactAngleFvPatchScalarField::
timeVaryingGammaContactAngleFvPatchScalarField
(
    const timeVaryingGammaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    gammaContactAngleFvPatchScalarField(gcpsf, iF),
    t0_(gcpsf.t0_),
    thetaT0_(gcpsf.thetaT0_),
    te_(gcpsf.te_),
    thetaTe_(gcpsf.thetaTe_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> timeVaryingGammaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    scalar t = patch().boundaryMesh().mesh().time().value();
    scalar theta0 = thetaT0_;

    if (t < t0_)
    {
        theta0 = thetaT0_;
    }
    else if (t > te_)
    {
        theta0 = thetaTe_;
    }
    else
    {
        theta0 = thetaT0_ + (t - t0_)*(thetaTe_ - thetaT0_)/(te_ - t0_);
    }

    return tmp<scalarField>(new scalarField(size(), theta0));
}


void timeVaryingGammaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("t0") << t0_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaT0") << thetaT0_ << token::END_STATEMENT << nl;
    os.writeKeyword("te") << te_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaTe") << thetaTe_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, timeVaryingGammaContactAngleFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
