/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Qro_(p.size(), 0.0),
    emissivity_(p.size(), 0.0)
{}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Qro_(ptf.Qro_),
    emissivity_(ptf.emissivity_)
{}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    Qro_("Qro", dict, p.size()),
    emissivity_("emissivity", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );

    }
    else
    {
         fvPatchScalarField::operator=(0.0);
    }
}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    Qro_(ptf.Qro_),
    emissivity_(ptf.emissivity_)
{}


Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    Qro_(ptf.Qro_),
    emissivity_(ptf.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
updateCoeffs()
{
    //Do nothing

    if (debug)
    {
        scalar Q = gSum((*this)*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << " heat[W]:" << Q
            << " wall radiative heat flux "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    Qro_.writeEntry("Qro", os);
    emissivity_.writeEntry("emissivity", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
    );
}
}


// ************************************************************************* //
