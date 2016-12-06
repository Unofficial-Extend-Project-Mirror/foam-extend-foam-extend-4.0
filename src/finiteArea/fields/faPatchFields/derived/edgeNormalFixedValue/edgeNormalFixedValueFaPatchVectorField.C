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

#include "edgeNormalFixedValueFaPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "areaFields.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(p, iF),
    refValue_(p.size(), 0)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& ptf,
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchVectorField(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchVectorField(p, iF, dict),
    refValue_("refValue", dict, p.size())
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& pivpvf
)
:
    fixedValueFaPatchVectorField(pivpvf),
    refValue_(pivpvf.refValue_)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& pivpvf,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(pivpvf, iF),
    refValue_(pivpvf.refValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeNormalFixedValueFaPatchVectorField::autoMap
(
    const faPatchFieldMapper& m
)
{
    fixedValueFaPatchVectorField::autoMap(m);
    refValue_.autoMap(m);
}


void Foam::edgeNormalFixedValueFaPatchVectorField::rmap
(
    const faPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFaPatchVectorField::rmap(ptf, addr);

    const edgeNormalFixedValueFaPatchVectorField& tiptf =
        refCast<const edgeNormalFixedValueFaPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::edgeNormalFixedValueFaPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Bug fix: update for moving mesh.  HJ, 15/Oct/2010
    operator==(refValue_*patch().edgeNormals());
}


void Foam::edgeNormalFixedValueFaPatchVectorField::write(Ostream& os) const
{
    fixedValueFaPatchVectorField::write(os);
    refValue_.writeEntry("refValue", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeFaPatchTypeField
(
    faPatchVectorField,
    edgeNormalFixedValueFaPatchVectorField
);

} // End namespace Foam


// ************************************************************************* //
