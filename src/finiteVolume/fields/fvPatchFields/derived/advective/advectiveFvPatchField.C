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

#include "advectiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "steadyStateDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
advectiveFvPatchField<Type>::advectiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    fieldInf_(pTraits<Type>::zero),
    lInf_(0.0),
    inletOutlet_(false),
    correctSupercritical_(false)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
advectiveFvPatchField<Type>::advectiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    fieldInf_(pTraits<Type>::zero),
    lInf_(0.0),
    inletOutlet_(dict.lookup("inletOutlet")),
    correctSupercritical_(dict.lookup("correctSupercritical"))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

    if (dict.readIfPresent("lInf", lInf_))
    {
        dict.lookup("fieldInf") >> fieldInf_;

        if (lInf_ < 0.0)
        {
            FatalIOErrorIn
            (
                "advectiveFvPatchField<Type>::"
                "advectiveFvPatchField"
                "(const fvPatch&, const Field<Type>&, const dictionary&)",
                dict
            )   << "unphysical lInf specified (lInf < 0)\n"
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalIOError);
        }
    }
}


template<class Type>
advectiveFvPatchField<Type>::advectiveFvPatchField
(
    const advectiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    fieldInf_(ptf.fieldInf_),
    lInf_(ptf.lInf_),
    inletOutlet_(ptf.inletOutlet_),
    correctSupercritical_(ptf.correctSupercritical_)
{}


template<class Type>
advectiveFvPatchField<Type>::advectiveFvPatchField
(
    const advectiveFvPatchField& ptpsf
)
:
    mixedFvPatchField<Type>(ptpsf),
    phiName_(ptpsf.phiName_),
    rhoName_(ptpsf.rhoName_),
    fieldInf_(ptpsf.fieldInf_),
    lInf_(ptpsf.lInf_),
    inletOutlet_(ptpsf.inletOutlet_),
    correctSupercritical_(ptpsf.correctSupercritical_)
{}


template<class Type>
advectiveFvPatchField<Type>::advectiveFvPatchField
(
    const advectiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptpsf, iF),
    phiName_(ptpsf.phiName_),
    rhoName_(ptpsf.rhoName_),
    fieldInf_(ptpsf.fieldInf_),
    lInf_(ptpsf.lInf_),
    inletOutlet_(ptpsf.inletOutlet_),
    correctSupercritical_(ptpsf.correctSupercritical_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<scalarField> advectiveFvPatchField<Type>::advectionSpeed() const
{
    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phip = this->patch().lookupPatchField
    (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop = this->patch().lookupPatchField
        (
            rhoName_,
            reinterpret_cast<const volScalarField*>(0),
            reinterpret_cast<const scalar*>(0)
        );

        return phip/(rhop*this->patch().magSf());
    }
    else
    {
        return phip/this->patch().magSf();
    }
}


template<class Type>
tmp<scalarField> advectiveFvPatchField<Type>::supercritical() const
{
    // In base class, the condition is never supercritical
    return tmp<scalarField>(new scalarField(this->size(), scalar(0)));
}


template<class Type>
void advectiveFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<Type, fvPatchField, volMesh>& field =
        this->db().objectRegistry::
        lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            this->dimensionedInternalField().name()
        );

    word ddtScheme
    (
        this->dimensionedInternalField().mesh().ddtScheme(field.name())
    );

    // Calculate the advection speed of the field wave
    // If the wave is incoming set the speed to 0.
    scalarField w = Foam::max(advectionSpeed(), scalar(0));

    // Get time-step
    scalarField deltaT(this->size(), this->db().time().deltaT().value());

    // For steady-state, calculate formal deltaT from the
    // under-relaxation factor
    if (ddtScheme == fv::steadyStateDdtScheme<Type>::typeName)
    {
        // Get under-relaxation factor
//         scalar urf =
//             this->dimensionedInternalField().mesh().relaxationFactor
//             (
//                 field.name()
//             );

//         deltaT = urf/(1 - urf)*1/(w*this->patch().deltaCoeffs() + SMALL);

        // Set delta t for Co = 1
        deltaT = 1/(w*this->patch().deltaCoeffs() + SMALL);
    }

    // Calculate the field wave coefficient alpha (See notes)
    scalarField alpha = w*deltaT*this->patch().deltaCoeffs();

    label patchi = this->patch().index();

    // Non-reflecting outflow boundary
    // If lInf_ defined setup relaxation to the value fieldInf_.
    // HJ, 25/Sep/2009
    if (lInf_ > SMALL)
    {
        // Calculate the field relaxation coefficient k (See notes)
        scalarField k = w*deltaT/lInf_;

        if
        (
            ddtScheme == fv::EulerDdtScheme<Type>::typeName
         || ddtScheme == fv::CrankNicholsonDdtScheme<Type>::typeName
        )
        {
            this->refValue() =
            (
                field.oldTime().boundaryField()[patchi] + k*fieldInf_
            )/(1.0 + k);

            this->valueFraction() = (1.0 + k)/(1.0 + alpha + k);
        }
        else if (ddtScheme == fv::backwardDdtScheme<Type>::typeName)
        {
            this->refValue() =
            (
                2.0*field.oldTime().boundaryField()[patchi]
              - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
              + k*fieldInf_
            )/(1.5 + k);

            this->valueFraction() = (1.5 + k)/(1.5 + alpha + k);
        }
        else if (ddtScheme == fv::steadyStateDdtScheme<Type>::typeName)
        {
            this->refValue() =
            (
                field.prevIter().boundaryField()[patchi] + k*fieldInf_
            )/(1.0 + k);

            this->valueFraction() = (1.0 + k)/(1.0 + alpha + k);
        }
        else
        {
            FatalErrorIn
            (
                "advectiveFvPatchField<Type>::updateCoeffs()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtScheme
                << "\n    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            ddtScheme == fv::EulerDdtScheme<Type>::typeName
         || ddtScheme == fv::CrankNicholsonDdtScheme<Type>::typeName
        )
        {
            this->refValue() = field.oldTime().boundaryField()[patchi];

            this->valueFraction() = 1.0/(1.0 + alpha);
        }
        else if (ddtScheme == fv::backwardDdtScheme<Type>::typeName)
        {
            this->refValue() =
            (
                2.0*field.oldTime().boundaryField()[patchi]
              - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
            )/1.5;

            this->valueFraction() = 1.5/(1.5 + alpha);
        }
        else if (ddtScheme == fv::steadyStateDdtScheme<Type>::typeName)
        {
            this->refValue() = field.prevIter().boundaryField()[patchi];

            this->valueFraction() = 1.0/(1.0 + alpha);
        }
        else
        {
            FatalErrorIn
            (
                "advectiveFvPatchField<Type>::updateCoeffs()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtScheme
                << "\n    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }

    // Get access to flux field
    fvsPatchField<scalar> phip = this->patch().lookupPatchField
    (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
    );

    // Treatment for supercritical inlet or outlet.  HJ, 28/Oct/2009
    // If the flow is faster than critical speed, the boundary condition
    // needs to behave as zero gradient at outflow and fixed value at inflow.
    // This is the case when eg.  a wave-transmissive pressure outlet
    // becomes supersonic.  Face selection for supercritical outlet is
    // done through a supercritical() virtual function: return 1 if
    // flow is supercritical and zero otherwise
    if (correctSupercritical_)
    {
        this->valueFraction() =
        (
            // If the flow is going out
            // - set zero gradient for supercritical; otherwise use current
            pos(phip)*
            (
                (scalar(1) - this->supercritical())
              + this->supercritical()*this->valueFraction()
            )
            // If the flow is going in
            // - set fixed value for supercritical; otherwise use current
          + neg(phip)*
            (
                this->supercritical()
              + (scalar(1) - this->supercritical())*this->valueFraction()
            )
        );

        this->refValue() =
            // If the flow is going out
            // - for supercritical, the value does not matter.  use fieldInf_
            // - for subcritical, use current value
            pos(phip)*
            (
                this->supercritical()*fieldInf_
              + (scalar(1) - this->supercritical())*this->refValue()
            )
            // If the flow is going in
            // - for supercritical use fieldInf_
          + neg(phip)*
            (
                this->supercritical()*fieldInf_
              + (scalar(1) - this->supercritical())*this->refValue()
            );
    }

    // Inlet-outlet treatment at the boundary.  HJ, 28/Oct/2009
    // If the flux is outgoing, the normal condition is used
    // if the flux is incoming, incoming value is set to fieldInf_
    // (value = fieldInf, valueFraction = 1)
    if (inletOutlet_)
    {
        this->refValue() = pos(phip)*this->refValue() + neg(phip)*fieldInf_;
        this->refGrad() = pTraits<Type>::zero;
        this->valueFraction() = pos(phip)*this->valueFraction() + neg(phip);
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void advectiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }

    if (lInf_ > SMALL)
    {
        os.writeKeyword("fieldInf") << fieldInf_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lInf") << lInf_
        << token::END_STATEMENT << endl;
    }

    os.writeKeyword("inletOutlet") << inletOutlet_
        << token::END_STATEMENT << nl;

    os.writeKeyword("correctSupercritical") << correctSupercritical_
            << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
