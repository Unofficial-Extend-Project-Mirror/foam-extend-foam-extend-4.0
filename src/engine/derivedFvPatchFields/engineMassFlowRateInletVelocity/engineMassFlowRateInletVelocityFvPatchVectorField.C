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

#include "engineMassFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "graph.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
engineMassFlowRateInletVelocityFvPatchVectorField::
engineMassFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}

Foam::
engineMassFlowRateInletVelocityFvPatchVectorField::
engineMassFlowRateInletVelocityFvPatchVectorField
(
    const engineMassFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}

Foam::
engineMassFlowRateInletVelocityFvPatchVectorField::
engineMassFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    phiName_("phi"),
    rhoName_("rho"),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand()),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{
    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    if (dict.found("rho"))
    {
        dict.lookup("rho") >> rhoName_;
    }
}

Foam::
engineMassFlowRateInletVelocityFvPatchVectorField::
engineMassFlowRateInletVelocityFvPatchVectorField
(
    const engineMassFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}

Foam::
engineMassFlowRateInletVelocityFvPatchVectorField::
engineMassFlowRateInletVelocityFvPatchVectorField
(
    const engineMassFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::engineMassFlowRateInletVelocityFvPatchVectorField::checkTable()
{
    if (!timeDataPtr_.valid())
    {
        timeDataPtr_.reset
        (
            new graph("title", "x", "y", IFstream(timeDataFileName_)())
        );
    }

    if (engineDB_.theta() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "engineTimeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()"
        )   << "current time (" << engineDB_.theta()
            << ") is less than the minimum in the data table ("
            << min(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the smallest time"
            << endl;
    }

    if (engineDB_.theta() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "engineTimeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()"
        )   << "current time (" << engineDB_.theta()
            << ") is greater than the maximum in the data table ("
            << max(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the largest time"
            << endl;
    }
}




void Foam::engineMassFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // a simpler way of doing this would be nice
    checkTable();
    
    scalar massFlowRate = 
    (
        interpolateXY
        (
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )
    );


    scalar avgU = -massFlowRate/gSum(patch().magSf());

    vectorField n = patch().nf();

    const surfaceScalarField& phi = db().lookupObject<surfaceScalarField>
    (
        phiName_
    );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // volumetric flow-rate
        operator==(n*avgU);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // mass flow-rate
        operator==(n*avgU/rhop);
    }
    else
    {
        FatalErrorIn
        (
            "engineMassFlowRateInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are incorrect"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::engineMassFlowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }

    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       engineMassFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
