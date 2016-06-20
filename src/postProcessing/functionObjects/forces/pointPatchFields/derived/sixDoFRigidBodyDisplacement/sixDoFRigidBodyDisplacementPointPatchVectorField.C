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

#include "sixDoFRigidBodyDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    motion_(),
    initialPoints_(p.localPoints()),
    rhoInf_(1.0)
{}


sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    motion_(dict),
    rhoInf_(readScalar(dict.lookup("rhoInf")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("initialPoints"))
    {
        initialPoints_ = vectorField("initialPoints", dict , p.size());
    }
    else
    {
        initialPoints_ = p.localPoints();
    }
}


sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const sixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_, mapper),
    rhoInf_(ptf.rhoInf_)
{}


sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const sixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_),
    rhoInf_(ptf.rhoInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sixDoFRigidBodyDisplacementPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);

    initialPoints_.autoMap(m);
}


void sixDoFRigidBodyDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const sixDoFRigidBodyDisplacementPointPatchVectorField& sDoFptf =
        refCast<const sixDoFRigidBodyDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchVectorField::rmap(sDoFptf, addr);

    initialPoints_.rmap(sDoFptf.initialPoints_, addr);
}


void sixDoFRigidBodyDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& ptPatch = this->patch();

    // Patch force data is valid for the current positions, so
    // calculate the forces on the motion object from this data, then
    // update the positions

    motion_.updatePosition(t.deltaTValue());

    dictionary forcesDict;

    forcesDict.add("patches", wordList(1, ptPatch.name()));
    forcesDict.add("rhoName", "rhoInf");
    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("CofR", motion_.centreOfMass());

    forces f("forces", db(), forcesDict);

    forces::forcesMoments fm = f.calcForcesMoment();

    // Get the forces on the patch faces at the current positions

    vector gravity = vector::zero;

    if (db().foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField g =
            db().lookupObject<uniformDimensionedVectorField>("g");

        gravity = g.value();
    }

    motion_.updateForce
    (
        fm.first().first() + fm.first().second() + gravity*motion_.mass(),
        fm.second().first() + fm.second().second(),
        t.deltaTValue()
    );

    Field<vector>::operator=
    (
        motion_.currentPosition(initialPoints_) - initialPoints_
    );

    fixedValuePointPatchVectorField::updateCoeffs();
}


void sixDoFRigidBodyDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    motion_.write(os);
    os.writeKeyword("rhoInf")
        << rhoInf_ << token::END_STATEMENT << nl;
    initialPoints_.writeEntry("initialPoints", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    sixDoFRigidBodyDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
