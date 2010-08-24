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

#include "fluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    reactivity_(0),
    gammaName_("gamma")
{
    this->gradient() = pTraits<Type>::zero;
}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    reactivity_(ptf.reactivity_),
    gammaName_(ptf.gammaName_)
{}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    reactivity_(readScalar(dict.lookup("reactivity"))),
    gammaName_(dict.lookup("gamma"))
{
    // Set dummy gradient
    this->gradient() = pTraits<Type>::zero;

    // Read the value entry from the dictionary
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "fluxFvPatchField<Type>::fluxFvPatchField"
            "("
            "const fvPatch& p,"
            "const DimensionedField<Type, volMesh>& iF,"
            "const dictionary& dict,"
            "const bool valueRequired"
            ")",
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    reactivity_(ptf.reactivity_),
    gammaName_(ptf.gammaName_)
{}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    reactivity_(ptf.reactivity_),
    gammaName_(ptf.gammaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Field<Type>& C = *this;

    const fvPatchField<scalar>& gammap =
        this->patch().lookupPatchField
        (
            gammaName_,
            reinterpret_cast<const volScalarField*>(NULL),
            reinterpret_cast<const scalar*>(NULL)
        );

    this->gradient() = reactivity_*C/gammap;

    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void fluxFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    os.writeKeyword("reactivity")
        << reactivity_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma")
        << gammaName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
