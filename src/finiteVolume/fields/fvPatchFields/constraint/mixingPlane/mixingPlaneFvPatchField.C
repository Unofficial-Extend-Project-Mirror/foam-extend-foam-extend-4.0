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

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlaneFvPatchField.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "volFields.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p))
{}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p))
{
    if (!isType<mixingPlaneFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not mixingPlane type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    if (!dict.found("value"))
    {
        // Grab the internal value for initialisation.
        fvPatchField<Type>::operator=(this->patchInternalField()());
    }

}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p))
{
    if (!isType<mixingPlaneFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField\n"
            "(\n"
            "    const mixingPlaneFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixingPlaneLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return neighbour field
template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::patchNeighbourField() const
{
    if(debug > 1)
    {
        Info << "mixingPlaneFvPatchField<Type>::patchNeighbourField(): for field: " << this->dimensionedInternalField().name();
        Info << " on patch: " << this->patch().name() << endl;
        Info << " surface Area: " << gSum(this->patch().magSf()) << endl;
    }

    const Field<Type>& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = mixingPlanePatch_.shadow().faceCells();

    Field<Type> sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = iField[sfc[i]];
    }

    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();

    if (this->size() > 0)
    {
        pnf = mixingPlanePatch_.interpolate(sField);
    }

    return tpnf;
}


template<class Type>
void mixingPlaneFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if(debug)
        Info << "Inside mixingPlaneFvPatchField<Type>::initEvaluate: for field: " << this->dimensionedInternalField().name() << endl;

    if(this->size() > 0)
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }

        Field<Type> pf
            (
                this->patch().weights()*this->patchInternalField()
                + (1.0 - this->patch().weights())*this->patchNeighbourField()
            );

        Field<Type>::operator=(pf);
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if(this->size() > 0)
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if(this->size() > 0)
    {
        // Get shadow face-cells and assemble shadow field
        const unallocLabelList& sfc = mixingPlanePatch_.shadow().faceCells();

        scalarField sField(sfc.size());

        forAll (sField, i)
        {
            sField[i] = psiInternal[sfc[i]];
        }

        scalarField pnf = mixingPlanePatch_.interpolate(sField);

        // Multiply the field by coefficients and add into the result
        const unallocLabelList& fc = mixingPlanePatch_.faceCells();

        forAll(fc, elemI)
        {
            result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{}

// Return averaged field on patch
template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::patchCircumferentialAverageField() const
{
    // Compute circum average of self
    return this->mixingPlanePatch_.circumferentialAverage(*this);
    
}

// Return the averaged field for patchInternalField
template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::patchInternalField() const
{
    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();

    if (this->size() > 0)
    {
        pnf = this->mixingPlanePatch_.circumferentialAverage(fvPatchField<Type>::patchInternalField());
    }

    return tpnf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
