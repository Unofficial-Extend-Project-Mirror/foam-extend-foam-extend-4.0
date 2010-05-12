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

#include "waveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("Undefined"),
    UName_("Undefined"),
    gamma_(0.0)
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.lookup("psi")),
    UName_(dict.lookup("U")),
    gamma_(readScalar(dict.lookup("gamma")))
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    UName_(ptf.UName_),
    gamma_(ptf.gamma_)
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    UName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    UName_(ptpsf.UName_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<scalarField> waveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip = this->patch().lookupPatchField
    (
        psiName_,
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
    );

    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>
        (this->phiName_);

    fvsPatchField<scalar> phip = this->patch().lookupPatchField
    (
        this->phiName_,
        reinterpret_cast<const surfaceScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
    );

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop = this->patch().lookupPatchField
        (
            this->rhoName_,
            reinterpret_cast<const volScalarField*>(NULL),
            reinterpret_cast<const scalar*>(NULL)
        );

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf() + sqrt(gamma_/psip);
}


template<class Type>
tmp<scalarField> waveTransmissiveFvPatchField<Type>::supercritical() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip = this->patch().lookupPatchField
    (
        psiName_,
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
    );

    const fvPatchVectorField& U =
        this->patch().lookupPatchField
        (
            UName_,
            reinterpret_cast<const volVectorField*>(NULL),
            reinterpret_cast<const vector*>(NULL)
        );

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return pos
    (
        mag(U.patchInternalField() & this->patch().Sf())/this->patch().magSf()
      - sqrt(gamma_/psip)
    );
}


template<class Type>
void waveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    advectiveFvPatchField<Type>::write(os);

    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
