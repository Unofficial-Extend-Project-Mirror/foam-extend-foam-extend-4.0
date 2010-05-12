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

#include "pressureDirectedInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    inletDir_(p.size())
{
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const pressureDirectedInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    inletDir_(ptf.inletDir_, mapper)
{}


pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    inletDir_("inletDirection", dict, p.size())
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


pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const pressureDirectedInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    inletDir_(pivpvf.inletDir_)
{}


pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const pressureDirectedInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureDirectedInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);

    inletDir_.autoMap(m);
}


void pressureDirectedInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const pressureDirectedInletOutletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureDirectedInletOutletVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void pressureDirectedInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    vectorField n = patch().nf();
    scalarField ndmagS = (n & inletDir_)*patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        refValue() = inletDir_*phip/ndmagS;
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>("rho");

        refValue() = inletDir_*phip/(rhop*ndmagS);
    }
    else
    {
        FatalErrorIn
        (
            "pressureDirectedInletOutletVelocityFvPatchVectorField::"
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
pressureDirectedInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureDirectedInletOutletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
