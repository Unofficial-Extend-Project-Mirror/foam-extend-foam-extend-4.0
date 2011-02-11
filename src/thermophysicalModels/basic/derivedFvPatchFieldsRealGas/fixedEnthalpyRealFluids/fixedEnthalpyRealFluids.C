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

#include "fixedEnthalpyRealFluids.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicPsiThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedEnthalpyRealFluids::fixedEnthalpyRealFluids
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::fixedEnthalpyRealFluids::fixedEnthalpyRealFluids
(
    const fixedEnthalpyRealFluids& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixedEnthalpyRealFluids::fixedEnthalpyRealFluids
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::fixedEnthalpyRealFluids::fixedEnthalpyRealFluids
(
    const fixedEnthalpyRealFluids& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::fixedEnthalpyRealFluids::fixedEnthalpyRealFluids
(
    const fixedEnthalpyRealFluids& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedEnthalpyRealFluids::updateCoeffs()
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

    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    Tw.evaluate();

    fvPatchScalarField& pw =
        const_cast<fvPatchScalarField&>(thermo.p().boundaryField()[patchi]);
    pw.evaluate();


    if (dimensionedInternalField().name() == "h")
    {
   
        operator==(thermo.h(pw,Tw, patchi));
    }
    else
    {
        operator==(thermo.hs(Tw, patchi));
    }

    fixedValueFvPatchScalarField::updateCoeffs();

    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedEnthalpyRealFluids
    );
}


// ************************************************************************* //
