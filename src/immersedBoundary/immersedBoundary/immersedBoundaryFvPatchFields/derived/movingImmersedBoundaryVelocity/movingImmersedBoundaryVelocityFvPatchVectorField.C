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

#include "movingImmersedBoundaryVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::updateIbValues()
{
    // Evaluate with complete boundary velocity
    vectorField::operator=
    (
        this->ibPatch().ibPolyPatch().motionDistance()/
        patch().boundaryMesh().mesh().time().deltaT().value()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    immersedBoundaryFieldBase<vector>(p, false, vector::zero)
{}


Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),   // Do not read data
    immersedBoundaryFieldBase<vector>
    (
        p,
        Switch(dict.lookup("setDeadValue")),
        vector(dict.lookup("deadValue"))
    )
{
    readPatchType(dict);
    updateIbValues();
}


Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const movingImmersedBoundaryVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF), // Do not map mixed data.  Set patchType later
    immersedBoundaryFieldBase<vector>(p, ptf.setDeadValue(), ptf.deadValue())
{
    // Copy the patch type since mixed data was not mapped
    this->setPatchType(ptf);

    updateIbValues();
}


Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const movingImmersedBoundaryVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    immersedBoundaryFieldBase<vector>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{}


Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const movingImmersedBoundaryVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    immersedBoundaryFieldBase<vector>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper&
)
{
    updateIbValues();
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField&,
    const labelList&
)
{
    if (size() != ibPatch().size())
    {
        updateIbValues();
    }
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::updateOnMotion()
{
    if (size() != ibPatch().size())
    {
        updateIbValues();
    }
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = dimensionedInternalField().mesh();

    if (mesh.changing())
    {
        const fvPatch& p = patch();

        // Get wall-parallel mesh motion velocity from immersed boundary
        vectorField Up = this->ibPatch().ibPolyPatch().motionDistance()/
            mesh.time().deltaT().value();

        const volVectorField& U =
            mesh.lookupObject<volVectorField>
            (
                dimensionedInternalField().name()
            );

        scalarField phip =
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

        // Warning: cannot use patch normal but the real face normal
        // THEY MAY NOT BE THE SAME!  HJ, 28/Mar/2019
        vectorField n = p.Sf()/(p.magSf());

        const scalarField& magSf = p.magSf();
        scalarField Un = phip/(magSf + VSMALL);

        // Adjust for surface-normal mesh motion flux
        vectorField::operator=(Up + n*(Un - (n & Up)));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    // Get non-constant reference to internal field
    vectorField& intField = const_cast<vectorField&>(this->internalField());

    // Set dead value
    this->setDeadValues(intField);

    // Evaluate mixed condition
    fixedValueFvPatchVectorField::evaluate();
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::manipulateMatrix
(
    fvVectorMatrix& matrix
)
{
    setDeadValues(matrix);
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);

    this->writeDeadData(os);

    vectorField::null().writeEntry("value", os);
    // writeEntry("value", os);

    this->writeField(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    movingImmersedBoundaryVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
