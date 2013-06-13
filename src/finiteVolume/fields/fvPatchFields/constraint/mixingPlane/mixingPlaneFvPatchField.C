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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void mixingPlaneFvPatchField<Type>::calcFluxMask() const
{
    // Find the flux field and calculate flux mask
    if (fluxAveraging_)
    {
        if (!this->db().objectRegistry::found(phiName_))
        {
            InfoIn
            (
                "void mixingPlaneFvPatchField<Type>::calcFluxMask()"
                ""
            )   << "Flux not found for flux averaging mixing plane on "
                << this->patch().name() << ".  Flux field: " << phiName_
                << endl;

            fluxMask_.setSize(mixingPlanePatch_.nProfileBands(), 0);
        }
        else
        {
            // Get flux field.  Used in decision for which faces
            // see the average and for flux-based weights
            const scalarField& phip = this->lookupPatchField
            (
                phiName_,
                reinterpret_cast<const surfaceScalarField*>(0),
                reinterpret_cast<const scalar*>(0)
            );

            if (mixingPlanePatch_.master())
            {
                // Calculate flux masks.  Master fluxes are used for the mask
                // Flux mask behaves like valueFraction:
                //     1 = use interpolated value
                //     0 = use patch internal field

                // Note: only master calculates the flux mask, as it is
                // contained on the profile.  Shadow will use the flux mask
                // from master.  HJ, 14/Feb/2013
                fluxMask_ = neg(mixingPlanePatch_.toProfile(phip));
            }
            else
            {
                // On the shadow side, flux mask is opposite from master
                fluxMask_ = 1 - shadowPatchField().fluxMask_;
            }

            // Calculate flux weights via circumferential average
            // Note: for zero flux, weights must be uniform and not zero
            // HJ, 4/Feb/2013

            // Note: flux weights exist on both master and slave
            scalar maxPhiP = 0;

            if (!phip.empty())
            {
                maxPhiP = max(phip);
            }

            reduce(maxPhiP, maxOp<scalar>());

            if (maxPhiP > SMALL)
            {
                fluxWeights_ = Foam::max(phip, scalar(0));
                fluxWeights_ /=
                    mixingPlanePatch_.circumferentialAverage(fluxWeights_)
                  + SMALL;
           }
            else
            {
                fluxWeights_ = 1;
            }
        }
    }
}


template<class Type>
const scalarField& mixingPlaneFvPatchField<Type>::fluxMask() const
{
    if (!this->updated())
    {
        // Force recalculation of flux masks
        calcFluxMask();        
    }

    return fluxMask_;
}


template<class Type>
const scalarField& mixingPlaneFvPatchField<Type>::fluxWeights() const
{
    if (!this->updated())
    {
        // Force recalculation of flux masks
        calcFluxMask();
    }

    return fluxWeights_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    fluxAveraging_(false),
    surfaceAveraging_(true),
    phiName_("phi"),
    fluxMask_(),
    fluxWeights_(p.size(), 0)
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
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    fluxAveraging_(dict.lookup("fluxAveraging")),
    surfaceAveraging_(dict.lookup("surfaceAveraging")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    fluxMask_(),
    fluxWeights_(p.size(), 0)
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
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    fluxAveraging_(ptf.fluxAveraging_),
    surfaceAveraging_(ptf.surfaceAveraging_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_, mapper)
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
    const mixingPlaneFvPatchField<Type>& ptf
)
:
    mixingPlaneLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch())),
    fluxAveraging_(ptf.fluxAveraging_),
    surfaceAveraging_(ptf.surfaceAveraging_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_)
{}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixingPlaneLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch())),
    fluxAveraging_(ptf.fluxAveraging_),
    surfaceAveraging_(ptf.surfaceAveraging_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return shadow field
template<class Type>
const mixingPlaneFvPatchField<Type>&
mixingPlaneFvPatchField<Type>::shadowPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->internalField()
        );

    return refCast<const mixingPlaneFvPatchField<Type> >
    (
        fld.boundaryField()[mixingPlanePatch_.shadowIndex()]
    );
}


template<class Type>
void mixingPlaneFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fluxMask_.setSize(0),
    fluxWeights_.autoMap(m);
}


template<class Type>
void mixingPlaneFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    fluxMask_.setSize(0);

    const mixingPlaneFvPatchField<Type>& tiptf =
        refCast<const mixingPlaneFvPatchField<Type> >(ptf);

    fluxWeights_.rmap(tiptf.fluxWeights_, addr);
}


template<class Type>
void mixingPlaneFvPatchField<Type>::updateCoeffs()
{
    // Force recalculation of flux masks
    calcFluxMask();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::patchNeighbourField() const
{
    // Get shadow patch internalField field
    Field<Type> sField = this->shadowPatchField().patchInternalField();

    if (fluxAveraging_)
    {
        // Flux averaging
        // - for outgoing flux, use zero gradient condition
        // - for incoming flux, use interpolated flux-weighted value

        const scalarField& mask = fluxMask();

        const scalarField& shadowFluxWeights =
            shadowPatchField().fluxWeights();

        // For outgoing flux, the value is identical to internal value
        // For incoming flux, calculate the average value of the
        // flux-weight shadow values coming out
        return
            mixingPlanePatch_.fromProfile
            (
                mask*mixingPlanePatch_.shadow().toProfile
                (
                    sField*shadowFluxWeights
                )
            )
          + mixingPlanePatch_.fromProfile(1 - mask)*this->patchInternalField();
    }
    else
    {
        // No flux averaging
        return mixingPlanePatch_.interpolate(sField);
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Evaluation is identical with or without flux averaging
    // as the masking and weighting is taken into account in
    // patchNeighbourField.  HJ, 4/Feb/2013

    const scalarField& w = this->patch().weights();            

    Field<Type>::operator=
    (
        w*this->patchInternalField()
      + (1 - w)*this->patchNeighbourField()
    );
}


template<class Type>
void mixingPlaneFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (fluxAveraging_)
    {
        // Flux averaging.
        // When flux averaging indicates outgoing flux,
        // use zero gradient condition: evaluating as internal field
        // When flux averaging indicates incoming flux,
        // use normal interpolation
        // Note that each face may have elements in multiple bands and
        // therefore interpolation from profile to patch is required
        // HJ, 11/Feb/2013
        const scalarField& mask = fluxMask();

        scalarField oneMFluxMask = mixingPlanePatch_.fromProfile(1 - mask);

        return pTraits<Type>::one*oneMFluxMask;
    }
    else
    {
        return Type(pTraits<Type>::one)*w;
    }
}

template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (fluxAveraging_)
    {
        // Flux averaging
        const scalarField& mask = fluxMask();

        scalarField fluxMask = mixingPlanePatch_.fromProfile(mask);

        return pTraits<Type>::one*fluxMask;
    }
    else
    {
        return Type(pTraits<Type>::one)*(1.0 - w);
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
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = mixingPlanePatch_.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    // Get local faceCells
    const unallocLabelList& fc = mixingPlanePatch_.faceCells();

    scalarField pnf = mixingPlanePatch_.interpolate(sField);

    // Multiply the field by coefficients and add into the result
    if (switchToLhs)
    {
        forAll(fc, elemI)
        {
            result[fc[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
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
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{}


template<class Type>
void mixingPlaneFvPatchField<Type>::patchInterpolate
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL
) const
{
    // Code moved from surfaceInterpolationScheme.C
    // HJ, 13/Jun/2013
    const label patchI = this->patch().index();

    // Use circumferential average of internal field when interpolating to
    // patch
    // HJ and MB, 13/Jun/2013
    if (surfaceAveraging_)
    {
        fField.boundaryField()[patchI] =
            pL*this->mixingPlanePatch_.circumferentialAverage
            (
                this->patchInternalField()
            )
          + (1 - pL)*this->patchNeighbourField();
    }
    else
    {
        fField.boundaryField()[patchI] =
            pL*this->patchInternalField()
          + (1 - pL)*this->patchNeighbourField();
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::patchInterpolate
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL,
    const scalarField& pY
) const
{
    // Code moved from surfaceInterpolationScheme.C
    // HJ, 13/Jun/2013
    const label patchI = this->patch().index();
    ::abort(); //HJ, HERE!!!
    // Use circumferential average of internal field
    // HJ and MB, 13/Jun/2013
    fField.boundaryField()[patchI] =
        pL*this->mixingPlanePatch_.circumferentialAverage
        (
            this->patchInternalField()
        )
      + pY*this->patchNeighbourField();
}


template<class Type>
void mixingPlaneFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeKeyword("fluxAveraging")
        << fluxAveraging_ << token::END_STATEMENT << nl;

    os.writeKeyword("surfaceAveraging")
        << surfaceAveraging_ << token::END_STATEMENT << nl;

    this->writeEntryIfDifferent(os, "phi", word("phi"), phiName_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
