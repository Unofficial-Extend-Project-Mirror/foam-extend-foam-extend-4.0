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

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Hrvoje Jasak, Wikki Ltd.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

\*---------------------------------------------------------------------------*/

#include "mixingPlaneEnthalpyJumpFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeTemplatePatchTypeField
    (
        fvPatchScalarField,
        mixingPlaneEnthalpyJumpFvPatchScalarField
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::mixingPlaneEnthalpyJumpFvPatchField<Foam::scalar>::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get access to relative and rotational velocity
    const word UrelName("Urel");
    const word UName("U");

    if
    (
        !this->db().objectRegistry::found(UrelName)
     || !this->db().objectRegistry::found(UName)
    )
    {
        // Velocities not available, do not update
        InfoIn
        (
            "void gradientEnthalpyFvPatchScalarField::"
            "updateCoeffs(const vectorField& Up)"
        )   << "Velocity fields " << UrelName << " or "
            << UName << " not found.  "
            << "Performing enthalpy value update" << endl;

        jump_ = 0;
    }
    else
    {
        const fvPatchVectorField& Urelp =
            lookupPatchField<volVectorField, vector>(UrelName);

        const fvPatchVectorField& Up =
            lookupPatchField<volVectorField, vector>(UName);

        if (rotating_)
        {
            jump_ =
                mag(Up.patchInternalField())*mag(Urelp.patchInternalField())
              - magSqr(Up.patchInternalField());
        }
        else
        {
            jump_ =
                mag(Up.patchNeighbourField())*mag(Urelp.patchNeighbourField())
              - magSqr(Up.patchNeighbourField());
        }
    }

    jumpMixingPlaneFvPatchField<scalar>::updateCoeffs();
}


// ************************************************************************* //
