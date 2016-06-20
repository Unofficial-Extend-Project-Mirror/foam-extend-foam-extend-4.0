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

#include "wedgeFvPatchVectorNFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeVectorTensorNWedgeFuncDefs(Type)                                  \
template<>                                                                    \
tmp<Field<Type> > wedgeFvPatchField<Type>::snGrad() const                     \
{                                                                             \
    return tmp<Field<Type> >                                                  \
    (                                                                         \
        new Field<Type>(size(), pTraits<Type>::zero)                          \
    );                                                                        \
}                                                                             \
                                                                              \
template<>                                                                    \
void wedgeFvPatchField<Type>::evaluate                                        \
(                                                                             \
    const Pstream::commsTypes commsType                                       \
)                                                                             \
{                                                                             \
    if (!updated())                                                           \
    {                                                                         \
        updateCoeffs();                                                       \
    }                                                                         \
                                                                              \
    this->operator==(patchInternalField());                                   \
}                                                                             \
                                                                              \
template<>                                                                    \
tmp<Field<Type> > wedgeFvPatchField<Type>::snGradTransformDiag()              \
const                                                                         \
{                                                                             \
    return tmp<Field<Type> >                                                  \
    (                                                                         \
        new Field<Type>(this->size(), pTraits<Type>::zero)                    \
    );                                                                        \
}


#define doMakePatchTypeField(type, Type, args...)                             \
                                                                              \
makeVectorTensorNWedgeFuncDefs(type)                                          \
                                                                              \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatch##Type##Field,                                                     \
    wedgeFvPatch##Type##Field                                                 \
);


forAllVectorNTypes(doMakePatchTypeField)

forAllTensorNTypes(doMakePatchTypeField)

forAllDiagTensorNTypes(doMakePatchTypeField)

forAllSphericalTensorNTypes(doMakePatchTypeField)


#undef doMakePatchTypeField

#undef makeVectorTensorNWedgeFuncDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
