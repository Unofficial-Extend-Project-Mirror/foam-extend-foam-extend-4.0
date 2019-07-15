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

#include "gradientEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
    const gradientEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
    const gradientEnthalpyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
    const gradientEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );

    const label patchi = patch().index();

    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);

    Tw.evaluate();

    if
    (
        dimensionedInternalField().name() == db().mangleFileName("h")
    )
    {
        gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
          + patch().deltaCoeffs()*
            (
                thermo.h(Tw, patchi)
              - thermo.h(Tw, patch().faceCells())
            );
    }
    else if
    (
        dimensionedInternalField().name() == db().mangleFileName("i")
    )
    {
        // Get access to relative and rotational velocity
        const word UName("U");
        const word URotName("URot");
        const word UThetaName("UTheta");

        if
        (
            !this->db().objectRegistry::foundObject<volVectorField>(UName)
         || !this->db().objectRegistry::foundObject<volVectorField>(URotName)
         || !this->db().objectRegistry::foundObject<volScalarField>(UThetaName)
        )
        {
             // Velocities not available, do not update
            InfoIn
            (
                "void gradientEnthalpyFvPatchScalarField::"
                "updateCoeffs(const vectorField& Up)"
            )   << "Velocity fields " << UName << " or "
                << URotName << " or "
                << UThetaName << " not found.  "
                << "Performing enthalpy value update for field "
                << this->dimensionedInternalField().name()
                << " and patch " << patchi
                << nl << "objects" << this->db().objectRegistry::sortedToc()
                << endl;

            gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
              + patch().deltaCoeffs()*
                (
                    thermo.h(Tw, patchi)
                  - thermo.h(Tw, patch().faceCells())
                );
        }
        else
        {
            const fvPatchVectorField& Up =
                lookupPatchField<volVectorField, vector>(UName);

            const fvPatchVectorField& URotp =
                lookupPatchField<volVectorField, vector>(URotName);

            const fvPatchScalarField& UThetap =
                lookupPatchField<volScalarField, scalar>(UThetaName);

            gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
              + patch().deltaCoeffs()*
                (
                    thermo.h(Tw, patchi)
                  - thermo.h(Tw, patch().faceCells())
                )
              + 0.5*mag(Up)*mag(Up.snGrad())
              - mag(UThetap)*mag(URotp.snGrad());
        }
    }
    else
    {
        gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
          + patch().deltaCoeffs()*
            (
                thermo.hs(Tw, patchi)
              - thermo.hs(Tw, patch().faceCells())
            );
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::gradientEnthalpyFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradientEnthalpyFvPatchScalarField
    );
}


// ************************************************************************* //
