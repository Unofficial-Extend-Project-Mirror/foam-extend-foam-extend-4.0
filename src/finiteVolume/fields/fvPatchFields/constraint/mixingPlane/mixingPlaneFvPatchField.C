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
void mixingPlaneFvPatchField<Type>::readMixingType() const
{
    const dictionary& dict =
        this->patch().boundaryMesh().mesh().schemesDict().subDict
        (
            "mixingPlane"
        );

    // Try reading field type
    word fieldName = this->dimensionedInternalField().name();

    if (dict.found(fieldName))
    {
        mixing_ = mixingPlaneInterpolation::mixingTypeNames_.read
        (
            dict.lookup(fieldName)
        );
    }
    else if (dict.found("default"))
    {
        mixing_ = mixingPlaneInterpolation::mixingTypeNames_.read
        (
            dict.lookup("default")
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::readMixingType() const",
            dict
        )   << "Cannot find mixing type for field "
            <<  this->dimensionedInternalField().name() << nl
            << "Please specify in fvSchemes in mixingPlane, "
            << "under field name " << nl
            << "Available types are "
            << mixingPlaneInterpolation::mixingTypeNames_.toc()
            << abort(FatalIOError);
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::calcFluxMask() const
{
    // Find the flux field and calculate flux mask
    if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    {
        if (!this->db().objectRegistry::found(phiName_))
        {
            InfoIn
            (
                "void mixingPlaneFvPatchField<Type>::calcFluxMask()"
                ""
            )   << "Flux not found for flux averaging mixing plane on "
                << this->patch().name() << " for field "
                << this->dimensionedInternalField().name()
                << ".  Flux field: " << phiName_
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
    mixing_(mixingPlaneInterpolation::MIXING_UNKNOWN),
    phiName_("phi"),
    fluxMask_(),
    fluxWeights_(p.size(), 0)
{
    // Cannot read mixing type here: internal field may not be set
    // HJ, 3/Jun/2015
}


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
    mixing_(mixingPlaneInterpolation::MIXING_UNKNOWN),
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

    // Read mixing type
    this->readMixingType();
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
    mixing_(ptf.mixing_),
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
    mixing_(ptf.mixing_),
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
    mixing_(ptf.mixing_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixingPlaneFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fluxMask_.clear();
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

    fluxMask_.clear();

    const mixingPlaneFvPatchField<Type>& tiptf =
        refCast<const mixingPlaneFvPatchField<Type> >(ptf);

    fluxWeights_.rmap(tiptf.fluxWeights_, addr);
}


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


// Return shadow field
template<class Type>
const mixingPlaneInterpolation::mixingType&
mixingPlaneFvPatchField<Type>::mixing() const
{
    // If mixing type is unknown, read it
    if (this->mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    return mixing_;
}


template<class Type>
void mixingPlaneFvPatchField<Type>::updateCoeffs()
{
    // Read mixing type
    readMixingType();

    // Force recalculation of flux masks
    calcFluxMask();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::patchNeighbourField() const
{
    // Get shadow patch internalField field
    Field<Type> sField = this->shadowPatchField().patchInternalField();

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        // Area-weighted averaging
        return mixingPlanePatch_.interpolate(sField);
    }
    else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
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
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        return this->patchInternalField();
    }
    else
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > mixingPlaneFvPatchField<Type>::"
            "patchNeighbourField() const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return this->patchInternalField();
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

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
        Field<Type>::operator=
        (
            w*this->mixingPlanePatch_.circumferentialAverage
            (
                this->patchInternalField()
            )
          + (1 - w)*this->patchNeighbourField()
        );
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        Field<Type>::operator=(this->patchInternalField());
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::initEvaluate\n"
            "(\n"
            "    const Pstream::commsTypes commsType\n"
            ")"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
    }
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
    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        return pTraits<Type>::one*w;
    }
    else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
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
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        return tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::one)
        );
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::valueInternalCoeffs\n"
            "(\n"
            "    const tmp<scalarField>& w\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        return pTraits<Type>::one*(1.0 - w);
    }
    else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    {
        // Flux averaging
        const scalarField& mask = fluxMask();

        scalarField fluxMask = mixingPlanePatch_.fromProfile(mask);

        return pTraits<Type>::one*fluxMask;
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        return  tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::zero)
        );
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::valueBoundaryCoeffs\n"
            "(\n"
            "    const tmp<scalarField>& w\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


// template<class Type>
// tmp<Field<Type> >
// mixingPlaneFvPatchField<Type>::gradientInternalCoeffs() const
// {
//     // If mixing type is unknown, read it
//     if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
//     {
//         this->readMixingType();
//     }
//
//     if
//     (
//         mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
//      || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
//     )
//     {
//         return -pTraits<Type>::one*this->patch().deltaCoeffs();
//     }
//     else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
//     {
//         return tmp<Field<Type> >
//         (
//             new Field<Type>(this->size(), pTraits<Type>::zero)
//         );
//     }
//     else
//     {
//         FatalErrorIn
//         (
//             "void mixingPlaneFvPatchField<Type>::gradientInternalCoeffs"
//             "gradientInternalCoeffs() const"
//         )   << "Unknown mixing type for patch " << this->patch().name()
//             << " for field "
//             <<  this->dimensionedInternalField().name()
//             << abort(FatalError);
//     }

//     // Dummy return to keep compiler happy
//     return tmp<Field<Type> >
//     (
//         new Field<Type>(this->size(), pTraits<Type>::zero)
//     );
// }


// template<class Type>
// tmp<Field<Type> >
// mixingPlaneFvPatchField<Type>::gradientBoundaryCoeffs() const
// {
//     // If mixing type is unknown, read it
//     if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
//     {
//         this->readMixingType();
//     }
//
//     if
//     (
//         mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
//      || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
//     )
//     {
//         -this->gradientInternalCoeffs();
//     }
//     else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
//     {
//         return tmp<Field<Type> >
//         (
//             new Field<Type>(this->size(), pTraits<Type>::zero)
//         );
//     }
//     else
//     {
//         FatalErrorIn
//         (
//             "void mixingPlaneFvPatchField<Type>::gradientBoundaryCoeffs"
//             "gradientBoundaryCoeffs() const"
//         )   << "Unknown mixing type for patch " << this->patch().name()
//             << " for field "
//             <<  this->dimensionedInternalField().name()
//             << abort(FatalError);
//     }

//     // Dummy return to keep compiler happy
//     return tmp<Field<Type> >
//     (
//         new Field<Type>(this->size(), pTraits<Type>::zero)
//     );
// }


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

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    // Use circumferential average of internal field when interpolating to
    // patch.  HJ and MB, 13/Jun/2013
    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
        fField.boundaryField()[patchI] =
            pL*this->mixingPlanePatch_.circumferentialAverage
            (
                this->patchInternalField()
            )
          + (1 - pL)*this->patchNeighbourField();
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        fField.boundaryField()[patchI] = this->patchInternalField();
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::patchInterpolate\n"
            "(\n"
            "    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,\n"
            "    const scalarField& pL\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
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

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    // Use circumferential average of internal field when interpolating to
    // patch.  HJ and MB, 13/Jun/2013
    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
        fField.boundaryField()[patchI] =
            pL*this->mixingPlanePatch_.circumferentialAverage
            (
                this->patchInternalField()
            )
          + pY*this->patchNeighbourField();
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        fField.boundaryField()[patchI] = this->patchInternalField();
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::patchInterpolate\n"
            "(\n"
            "    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,\n"
            "    const scalarField& pL,\n"
            "    const scalarField& pY\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::patchFlux
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<Type>& matrix
) const
{
    const label patchI = this->patch().index();

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
//         flux.boundaryField()[patchI] =
//             cmptMultiply
//             (
//                 matrix.internalCoeffs()[patchI],
//                 this->mixingPlanePatch_.circumferentialAverage
//                 (
//                     this->patchInternalField()
//                 )
//             )
//           - cmptMultiply
//             (
//                 matrix.boundaryCoeffs()[patchI],
//                 mixingPlanePatch_.shadow().circumferentialAverage
//                 (
//                     this->patchNeighbourField()
//                 )
//             );

        //  Alternative, not sure which is correct
        flux.boundaryField()[patchI] =
            cmptMultiply
            (
                matrix.internalCoeffs()[patchI],
                this->patchInternalField()
            )
          - cmptMultiply
            (
                matrix.boundaryCoeffs()[patchI],
                this->patchNeighbourField()
            );
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        flux.boundaryField()[patchI] = pTraits<Type>::zero;
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::patchFlux\n"
            "(\n"
            "    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,\n"
            "    const fvMatrix<Type>& matrix\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
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

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
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
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        // Do nothing
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::initInterfaceMatrixUpdate\n"
            "(\n"
            "    const scalarField& psiInternal,\n"
            "    scalarField& result,\n"
            "    const lduMatrix&,\n"
            "    const scalarField& coeffs,\n"
            "    const direction cmpt,\n"
            "    const Pstream::commsTypes commsType,\n"
            "    const bool switchToLhs\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->dimensionedInternalField().name()
            << abort(FatalError);
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
void mixingPlaneFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    this->writeEntryIfDifferent(os, "phi", word("phi"), phiName_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
