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

#include "mixedInternalEnergyRealFluids.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicPsiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixedInternalEnergyRealFluids::mixedInternalEnergyRealFluids
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


mixedInternalEnergyRealFluids::mixedInternalEnergyRealFluids
(
    const mixedInternalEnergyRealFluids& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


mixedInternalEnergyRealFluids::mixedInternalEnergyRealFluids
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


mixedInternalEnergyRealFluids::mixedInternalEnergyRealFluids
(
    const mixedInternalEnergyRealFluids& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


mixedInternalEnergyRealFluids::mixedInternalEnergyRealFluids
(
    const mixedInternalEnergyRealFluids& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mixedInternalEnergyRealFluids::updateCoeffs()
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
    refValue() = thermo.e(pw,Tw.refValue(), patchi);
    refGrad() = thermo.Cv(pw,Tw, patchi)*Tw.refGrad()
      + patch().deltaCoeffs()*
        (
            thermo.e(pw,Tw, patchi)
          - thermo.e(pw,Tw, patch().faceCells())
        );

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mixedInternalEnergyRealFluids);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
