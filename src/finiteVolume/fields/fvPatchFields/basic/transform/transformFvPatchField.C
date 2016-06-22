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

#include "transformFvPatchField.H"
#include "IOstreams.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
transformFvPatchField<Type>::transformFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
   fvPatchField<Type>(p, iF)
{}


template<class Type>
transformFvPatchField<Type>::transformFvPatchField
(
    const transformFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
transformFvPatchField<Type>::transformFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict)
{}


template<class Type>
transformFvPatchField<Type>::transformFvPatchField
(
    const transformFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf)
{}


template<class Type>
transformFvPatchField<Type>::transformFvPatchField
(
    const transformFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > transformFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return pTraits<Type>::one - this->snGradTransformDiag();
}


template<class Type>
tmp<Field<Type> > transformFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        *this
      - cmptMultiply
        (
            this->valueInternalCoeffs(this->patch().weights()),
            this->patchInternalField()
        );
}


template<class Type>
tmp<Field<Type> > transformFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -this->patch().deltaCoeffs()*this->snGradTransformDiag();
}


template<class Type>
tmp<Field<Type> > transformFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        this->snGrad()
      - cmptMultiply
        (
            this->gradientInternalCoeffs(),
            this->patchInternalField()
        );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void transformFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    this->evaluate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
