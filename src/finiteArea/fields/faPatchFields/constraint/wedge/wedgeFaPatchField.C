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

#include "wedgeFaPatch.H"
#include "wedgeFaPatchField.H"
#include "transformField.H"
#include "symmTransform.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFaPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "wedgeFaPatchField<Type>::wedgeFaPatchField\n"
            "(\n"
            "    const wedgeFaPatchField<Type>& ptf,\n"
            "    const faPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const faPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    if (!isType<wedgeFaPatch>(p))
    {
        FatalIOErrorIn
        (
            "wedgeFaPatchField<Type>::wedgeFaPatchField\n"
            "(\n"
            "    const faPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not wedge type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    evaluate();
}


template<class Type>
wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return gradient at boundary
template<class Type>
tmp<Field<Type> > wedgeFaPatchField<Type>::snGrad() const
{
    Field<Type> pif = this->patchInternalField();
    return
    (
        transform(refCast<const wedgeFaPatch>(this->patch()).faceT(), pif)
      - pif
    )*(0.5*this->patch().deltaCoeffs());
}


// Evaluate the patch field
template<class Type>
void wedgeFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    faPatchField<Type>::operator==
    (
        transform
        (
            refCast<const wedgeFaPatch>(this->patch()).edgeT(),
            this->patchInternalField()
        )
    );
}


// Return defining fields
template<class Type>
tmp<Field<Type> > wedgeFaPatchField<Type>::snGradTransformDiag() const
{
    const diagTensor diagT =
        0.5*diag(I - refCast<const wedgeFaPatch>(this->patch()).faceT());

    const vector diagV(diagT.xx(), diagT.yy(), diagT.zz());

    return tmp<Field<Type> >
    (
        new Field<Type>
        (
            this->size(),
            transformMask<Type>
            (
                pow
                (
                    diagV,
                    pTraits<typename powProduct<vector, pTraits<Type>::rank>
                    ::type>::zero
                )
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
