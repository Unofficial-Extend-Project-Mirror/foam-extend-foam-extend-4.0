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

#include "mixedEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedEnthalpyFvPatchScalarField::mixedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixedEnthalpyFvPatchScalarField::mixedEnthalpyFvPatchScalarField
(
    const mixedEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedEnthalpyFvPatchScalarField::mixedEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedEnthalpyFvPatchScalarField::mixedEnthalpyFvPatchScalarField
(
    const mixedEnthalpyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixedEnthalpyFvPatchScalarField::mixedEnthalpyFvPatchScalarField
(
    const mixedEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedEnthalpyFvPatchScalarField::updateCoeffs()
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

    mixedFvPatchScalarField& Tw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi])
    );

    Tw.evaluate();

    valueFraction() = Tw.valueFraction();

    if
    (
        dimensionedInternalField().name() == db().mangleFileName("h")
    )
    {
        refValue() = thermo.h(Tw.refValue(), patchi);
        refGrad() = thermo.Cp(Tw, patchi)*Tw.refGrad()
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
            !this->db().objectRegistry::found(UName)
         || !this->db().objectRegistry::found(URotName)
         || !this->db().objectRegistry::found(UThetaName)
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
                << endl;

             refValue() = thermo.h(Tw.refValue(), patchi);
             refGrad() = thermo.Cp(Tw, patchi)*Tw.refGrad()
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

            refValue() = thermo.h(Tw.refValue(), patchi)
                + 0.5*magSqr(Up)
                - mag(UThetap)*mag(URotp);

            refGrad() = thermo.Cp(Tw, patchi)*Tw.refGrad()
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
        refValue() = thermo.hs(Tw.refValue(), patchi);
        refGrad() = thermo.Cp(Tw, patchi)*Tw.refGrad()
          + patch().deltaCoeffs()*
            (
                thermo.hs(Tw, patchi)
              - thermo.hs(Tw, patch().faceCells())
            );
    }

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedEnthalpyFvPatchScalarField
    );
}


// ************************************************************************* //
