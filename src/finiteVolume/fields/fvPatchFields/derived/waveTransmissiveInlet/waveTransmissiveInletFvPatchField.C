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

#include "waveTransmissiveInletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
waveTransmissiveInletFvPatchField<Type>::waveTransmissiveInletFvPatchField
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
waveTransmissiveInletFvPatchField<Type>::waveTransmissiveInletFvPatchField
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
waveTransmissiveInletFvPatchField<Type>::waveTransmissiveInletFvPatchField
(
    const waveTransmissiveInletFvPatchField& ptf,
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
waveTransmissiveInletFvPatchField<Type>::waveTransmissiveInletFvPatchField
(
    const waveTransmissiveInletFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    UName_(ptpsf.UName_),
    gamma_(ptpsf.gamma_)
{}


template<class Type>
waveTransmissiveInletFvPatchField<Type>::waveTransmissiveInletFvPatchField
(
    const waveTransmissiveInletFvPatchField& ptpsf,
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
tmp<scalarField>
 waveTransmissiveInletFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip = this->patch().lookupPatchField
    (
        psiName_,
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
    );

    const surfaceScalarField& phi =
        this->db().objectRegistry::template lookupObject<surfaceScalarField>
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
    // by subtracting the speed of sound (sqrt(gamma_/psi)) from the
    // component of the velocity normal to the boundary
    // Note: at the inlet, the flux has got a negative sign
    // Bug fix: wrong speed.  HJ, 26/Apr/2010
    return sqrt(gamma_/psip) + phip/this->patch().magSf();
}


template<class Type>
tmp<scalarField>
waveTransmissiveInletFvPatchField<Type>::supercritical() const
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
void waveTransmissiveInletFvPatchField<Type>::write(Ostream& os) const
{
    advectiveFvPatchField<Type>::write(os);
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
