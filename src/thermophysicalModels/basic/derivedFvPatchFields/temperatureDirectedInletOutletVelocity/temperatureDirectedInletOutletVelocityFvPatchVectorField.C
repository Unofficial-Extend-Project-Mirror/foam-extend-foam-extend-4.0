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

#include "temperatureDirectedInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureDirectedInletOutletVelocityFvPatchVectorField::
temperatureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    TName_("T"),
    T0_(p.size(), 0.0),
    inletDir_(p.size()),
    cylindricalCCS_(0),
    omega_(vector::zero)
{
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


temperatureDirectedInletOutletVelocityFvPatchVectorField::
temperatureDirectedInletOutletVelocityFvPatchVectorField
(
    const temperatureDirectedInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    TName_(ptf.TName_),
    T0_(ptf.T0_, mapper),
    inletDir_(ptf.inletDir_, mapper),
    cylindricalCCS_(ptf.cylindricalCCS_),
    omega_(ptf.omega_)
{}


temperatureDirectedInletOutletVelocityFvPatchVectorField::
temperatureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    TName_("T"),
    T0_("T0", dict, p.size()),
    inletDir_("inletDirection", dict, p.size()),
    cylindricalCCS_(dict.lookup("cylindricalCCS")),
    omega_(dict.lookup("omega"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


temperatureDirectedInletOutletVelocityFvPatchVectorField::
temperatureDirectedInletOutletVelocityFvPatchVectorField
(
    const temperatureDirectedInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    TName_(pivpvf.TName_),
    T0_(pivpvf.T0_),
    inletDir_(pivpvf.inletDir_),
    cylindricalCCS_(pivpvf.cylindricalCCS_),
    omega_(pivpvf.omega_)
{}


temperatureDirectedInletOutletVelocityFvPatchVectorField::
temperatureDirectedInletOutletVelocityFvPatchVectorField
(
    const temperatureDirectedInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    TName_(pivpvf.TName_),
    T0_(pivpvf.T0_),
    inletDir_(pivpvf.inletDir_),
    cylindricalCCS_(pivpvf.cylindricalCCS_),
    omega_(pivpvf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureDirectedInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);

    inletDir_.autoMap(m);
}


void temperatureDirectedInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const temperatureDirectedInletOutletVelocityFvPatchVectorField& tiptf =
        refCast<const temperatureDirectedInletOutletVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void temperatureDirectedInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& C = patch().Cf();
    vectorField rotationVelocity = omega_ ^ C;
    vectorField inletDirCCS = inletDir_/mag(inletDir_);
    vectorField inletDir(inletDirCCS);

    if (cylindricalCCS_)
    {
        scalar radius;
        forAll(C, facei)
        {
            radius = sqrt(sqr(C[facei].x())+ sqr(C[facei].y()));
            inletDir[facei].z() = inletDirCCS[facei].z();

            if (radius > 0.0)
            {
                inletDir[facei].x() =
                    (C[facei].x()*inletDirCCS[facei].x()
                    - C[facei].y()*inletDirCCS[facei].y())/radius;

                inletDir[facei].y() =
                    (C[facei].y()*inletDirCCS[facei].x()
                    + C[facei].x()*inletDirCCS[facei].y())/radius;
            }
            else
            {
                inletDir[facei].x() = 0.0;
                inletDir[facei].y() = 0.0;
            }
        }
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        vectorField n = patch().nf();
        scalarField ndmagS = (n & inletDir_)*patch().magSf();

        refValue() = inletDir*phip/ndmagS - rotationVelocity;
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        volScalarField Cp = thermo.Cp();

        const fvPatchField<scalar>& Cpp =
            patch().patchField<volScalarField, scalar>(Cp);

        refValue() =
            inletDir*sqrt(2.0*Cpp*max(T0_-Tp,VSMALL)) - rotationVelocity;
    }
    else
    {
        FatalErrorIn
        (
            "temperatureDirectedInletOutletVelocityFvPatchVectorField::"
            "updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    valueFraction() = 1.0 - pos(phip);

    mixedFvPatchVectorField::updateCoeffs();
}


void
temperatureDirectedInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    T0_.writeEntry("T0", os);
    inletDir_.writeEntry("inletDirection", os);
    os.writeKeyword("cylindricalCCS") << cylindricalCCS_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega")<< omega_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    temperatureDirectedInletOutletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
