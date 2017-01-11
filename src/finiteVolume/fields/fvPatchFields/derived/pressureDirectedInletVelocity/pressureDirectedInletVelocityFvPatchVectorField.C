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

#include "pressureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    inletDir_(p.size(), vector(1, 0, 0))
{}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletDir_(ptf.inletDir_, mapper)
{}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (!inletDir_.empty())
    {
        if (min(mag(inletDir_)) < SMALL)
        {
            FatalErrorIn
            (
                "pressureDirectedInletVelocityFvPatchVectorField::\n"
                "pressureDirectedInletVelocityFvPatchVectorField\n"
                "(\n"
                "    const fvPatch& p,\n"
                "    const DimensionedField<vector, volMesh>& iF,\n"
                "    const dictionary& dict\n"
                ")"
            )    << "Badly defined inlet direction for field "
                 << this->dimensionedInternalField().name()
                 << " and patch " << this->patch().name()
                 << abort(FatalError);
        }
    }       

    // Normalise to obtain the flow direction
    inletDir_ /= (mag(inletDir_) + SMALL);
}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_)
{}


pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureDirectedInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
}


void pressureDirectedInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const pressureDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureDirectedInletVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void pressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!this->db().objectRegistry::found(phiName_))
    {
        // Flux not available, do not update
        InfoIn
        (
            "void pressureDirectedInletVelocityFvPatchVectorField::"
            "updateCoeffs()"
        )   << "Flux field " << phiName_ << " not found.  "
            << "Performing fixed value update" << endl;

        fixedValueFvPatchVectorField::updateCoeffs();

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
        operator==(inletDir_*phip/ndmagS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        if (!this->db().objectRegistry::found(rhoName_))
        {
            // Rho not available, do not update
            fixedValueFvPatchVectorField::updateCoeffs();

            return;
        }

        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(inletDir_*phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorIn
        (
            "pressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void pressureDirectedInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void pressureDirectedInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureDirectedInletVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
