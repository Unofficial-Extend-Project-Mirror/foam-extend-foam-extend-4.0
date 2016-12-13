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

#include "IOobject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(ptf, mapper),
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else if (!valueRequired)
    {
        fvPatchField<Type>::operator=(pTraits<Type>::zero);
    }
    else
    {
        FatalIOErrorIn
        (
            "fvPatchField<Type>::fvPatchField"
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
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false),
    patchType_(ptf.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::fvPatchField<Type>::db() const
{
    //HR 12.3.10: Lookup fields from the field DB rather than the mesh
    return internalField_.db();
}


template<class Type>
template<class GeometricField, class Type2>
const typename GeometricField::PatchFieldType&
Foam::fvPatchField<Type>::lookupPatchField
(
    const word& name,
    const GeometricField*,
    const Type2*
) const
{
    return patch_.patchField<GeometricField, Type2>
    (
        internalField_.db().objectRegistry::template lookupObject<GeometricField>(name)
    );
}


template<class Type>
void Foam::fvPatchField<Type>::check(const fvPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorIn("PatchField<Type>::check(const fvPatchField<Type>&)")
            << "different patches for fvPatchField<Type>s"
            << abort(FatalError);
    }
}


// Return gradient at boundary
template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fvPatchField<Type>::snGrad() const
{
    return (*this - patchInternalField())*patch_.deltaCoeffs();
}


// Return internal field next to patch as patch field
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::fvPatchField<Type>::patchInternalField() const
{
    return patch_.patchInternalField(internalField_);
}


template<class Type>
void Foam::fvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


template<class Type>
void Foam::fvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated_)
    {
        updateCoeffs();
    }

    updated_ = false;
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix(fvMatrix<Type>& matrix)
{
    // do nothing
}


template<class Type>
void Foam::fvPatchField<Type>::patchInterpolate
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL
) const
{
    // Code moved from surfaceInterpolationScheme.C
    // HJ, 29/May/2013
    const label patchI = this->patch().index();

     // Virtual function for patch face interpolate.  HJ, 13/Jun/2013
     if (this->coupled())
     {
         // Coupled patch
         fField.boundaryField()[patchI] =
             pL*this->patchInternalField()
           + (1 - pL)*this->patchNeighbourField();
     }
     else
     {
         // Uncoupled patch, re-use face values
         fField.boundaryField()[patchI] = *this;
     }
}


template<class Type>
void Foam::fvPatchField<Type>::patchInterpolate
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL,
    const scalarField& pY
) const
{
    // Code moved from surfaceInterpolationScheme.C
    // HJ, 29/May/2013
    const label patchI = this->patch().index();

    // Virtual function for patch face interpolate.  HJ, 13/Jun/2013
    if (this->coupled())
    {
        // Coupled patch
        fField.boundaryField()[patchI] =
            pL*this->patchInternalField()
          + pY*this->patchNeighbourField();
    }
    else
    {
        // Uncoupled patch, re-used face values
        fField.boundaryField()[patchI] = *this;
    }
}


template<class Type>
void Foam::fvPatchField<Type>::patchFlux
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<Type>& matrix
) const
{
    // Code moved from fvMatrix.C
    // HJ, 29/May/2013
    const label patchI = this->patch().index();

    // Virtual function for patch flux.  HJ, 29/May/2013
    if (this->coupled())
    {
        // Coupled patch
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
    else
    {
        // Uncoupled patch
        flux.boundaryField()[patchI] =
            cmptMultiply
            (
                matrix.internalCoeffs()[patchI],
                this->patchInternalField()
            )
          - matrix.boundaryCoeffs()[patchI];
    }
}


template<class Type>
void Foam::fvPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;

    if (patchType_.size())
    {
        os.writeKeyword("patchType") << patchType_
            << token::END_STATEMENT << nl;
    }
}


template<class Type>
template<class EntryType>
void Foam::fvPatchField<Type>::writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const EntryType& value1,
    const EntryType& value2
) const
{
    if (value1 != value2)
    {
        os.writeKeyword(entryName) << value2 << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator*=(const fvPatchField<scalar>& ptf)"
        )   << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator/=(const fvPatchField<scalar>& ptf)"
        )   << "    incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const fvPatchField<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newFvPatchField.C"

// ************************************************************************* //
