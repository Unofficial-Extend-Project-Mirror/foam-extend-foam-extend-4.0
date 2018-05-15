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

#include "movingImmersedBoundaryVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::setDeadValues()
{
    // Fix the value in dead cells
    if (setDeadValue_)
    {
        const labelList& dc = ibPatch_.ibPolyPatch().deadCells();

        // Get non-const access to internal field
        vectorField& psiI = const_cast<vectorField&>(this->internalField());

        forAll (dc, dcI)
        {
            psiI[dc[dcI]] = deadValue_;
        }
    }
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
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    setDeadValue_(false),
    deadValue_(vector::zero)
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
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    setDeadValue_(dict.lookup("setDeadValue")),
    deadValue_(dict.lookup("deadValue"))
{
    // Evaluate with complete boundary velocity
    vectorField::operator=
    (
        ibPatch_.ibPolyPatch().motionDistance()/
        patch().boundaryMesh().mesh().time().deltaT().value()
    );
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
    fixedValueFvPatchVectorField(p, iF),   // Do not  data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{
    // Evaluate with complete boundary velocity
    vectorField::operator=
    (
        ibPatch_.ibPolyPatch().motionDistance()/
        patch().boundaryMesh().mesh().time().deltaT().value()
    );
}


Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const movingImmersedBoundaryVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    ibPatch_(ptf.ibPatch()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{}


Foam::movingImmersedBoundaryVelocityFvPatchVectorField::
movingImmersedBoundaryVelocityFvPatchVectorField
(
    const movingImmersedBoundaryVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    ibPatch_(ptf.ibPatch()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
        vectorField Up = ibPatch_.ibPolyPatch().motionDistance()/
            mesh.time().deltaT().value();

        const volVectorField& U =
            mesh.lookupObject<volVectorField>
            (
                dimensionedInternalField().name()
            );

        scalarField phip =
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

        vectorField n = p.nf();
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
    // Set dead value
    this->setDeadValues();

    // Evaluate mixed condition
    fixedValueFvPatchVectorField::evaluate();
}


void Foam::movingImmersedBoundaryVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("patchType")
        << immersedBoundaryFvPatch::typeName << token::END_STATEMENT << nl;
    os.writeKeyword("setDeadValue")
        << setDeadValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("deadValue")
        << deadValue_ << token::END_STATEMENT << nl;

    vectorField::null().writeEntry("value", os);
    // writeEntry("value", os);

    // Write VTK on master only
    if (Pstream::master())
    {
        // Add parallel reduction of all faces and data to proc 0
        // and write the whola patch together

        // Write immersed boundary data as a vtk file
        autoPtr<surfaceWriter<vector> > writerPtr =
            surfaceWriter<vector>::New("vtk");

        // Get the intersected patch
        const standAlonePatch& ts = ibPatch_.ibPolyPatch().ibPatch();

        writerPtr->write
        (
            this->dimensionedInternalField().path(),
            ibPatch_.name(),
            ts.points(),
            ts,
            this->dimensionedInternalField().name(),
            *this,
            surfaceWriterBase::FACE_DATA
        );
    }
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
