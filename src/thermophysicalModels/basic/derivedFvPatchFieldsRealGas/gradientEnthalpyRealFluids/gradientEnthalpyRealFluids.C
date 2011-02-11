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

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "gradientEnthalpyRealFluids.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicPsiThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientEnthalpyRealFluids::gradientEnthalpyRealFluids
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::gradientEnthalpyRealFluids::gradientEnthalpyRealFluids
(
    const gradientEnthalpyRealFluids& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::gradientEnthalpyRealFluids::gradientEnthalpyRealFluids
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::gradientEnthalpyRealFluids::gradientEnthalpyRealFluids
(
    const gradientEnthalpyRealFluids& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::gradientEnthalpyRealFluids::gradientEnthalpyRealFluids
(
    const gradientEnthalpyRealFluids& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientEnthalpyRealFluids::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicPsiThermo& thermo = db().lookupObject<basicPsiThermo>
    (
        "thermophysicalProperties"
    );

    const label patchi = patch().index();

    fvPatchScalarField& pw =
        const_cast<fvPatchScalarField&>(thermo.p().boundaryField()[patchi]);

    pw.evaluate();

    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);

    Tw.evaluate();

    if (dimensionedInternalField().name() == "h")
    {
        gradient() = thermo.Cp(pw,Tw, patchi)*Tw.snGrad()
        + patch().deltaCoeffs()*
        (
            thermo.h(pw,Tw, patchi)
          - thermo.h(pw,Tw, patch().faceCells())
        );
    }
 /*   else
    {
        gradient() = thermo.Cp(Tw, patchi)*Tw.snGrad()
        + patch().deltaCoeffs()*
        (
            thermo.hs(Tw, patchi)
          - thermo.hs(Tw, patch().faceCells())
        );
    }
*/
    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradientEnthalpyRealFluids
    );
}


// ************************************************************************* //
