/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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

#include "mixedEnthalpyRealFluids.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicPsiThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedEnthalpyRealFluids::mixedEnthalpyRealFluids
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


Foam::mixedEnthalpyRealFluids::mixedEnthalpyRealFluids
(
    const mixedEnthalpyRealFluids& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedEnthalpyRealFluids::mixedEnthalpyRealFluids
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedEnthalpyRealFluids::mixedEnthalpyRealFluids
(
    const mixedEnthalpyRealFluids& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixedEnthalpyRealFluids::mixedEnthalpyRealFluids
(
    const mixedEnthalpyRealFluids& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedEnthalpyRealFluids::updateCoeffs()
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

    mixedFvPatchScalarField& Tw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi])
    );

    Tw.evaluate();

    fvPatchScalarField& pw =
        const_cast<fvPatchScalarField&>(thermo.p().boundaryField()[patchi]);
    pw.evaluate();

    valueFraction() = Tw.valueFraction();

    if (dimensionedInternalField().name() == "h")
    {
        refValue() = thermo.h(pw,Tw.refValue(), patchi);
        refGrad() = thermo.Cp(pw,Tw, patchi)*Tw.refGrad()
        + patch().deltaCoeffs()*
         (
            thermo.h(pw,Tw, patchi)
          - thermo.h(pw,Tw, patch().faceCells())
         );
    }
  /*  else
    {
        refValue() = thermo.hs(Tw.refValue(), patchi);
        refGrad() = thermo.CpBC(pw,Tw, patchi)*Tw.refGrad()
        + patch().deltaCoeffs()*
         (
            thermo.hs(Tw, patchi)
          - thermo.hs(Tw, patch().faceCells())
         );
    }*/

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedEnthalpyRealFluids
    );
}


// ************************************************************************* //
