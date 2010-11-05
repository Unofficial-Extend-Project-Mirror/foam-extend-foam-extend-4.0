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

#include "wedgeFvPatchVectorNFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeVectorTensorNWedgeFuncDefs(Type)                    \
template<>                                                      \
tmp<Field<Type> > wedgeFvPatchField<Type>::snGrad() const       \
{                                                               \
    return tmp<Field<Type> >                                    \
    (                                                           \
        new Field<Type>(size(), pTraits<Type>::zero)            \
    );                                                          \
}                                                               \
                                                                \
template<>                                                      \
void wedgeFvPatchField<Type>::evaluate(                         \
    const Pstream::commsTypes commsType                         \
)                                                               \
{                                                               \
    if (!updated())                                             \
    {                                                           \
        updateCoeffs();                                         \
    }                                                           \
                                                                \
    operator==(patchInternalField());                           \
}                                                               \
                                                                \
template<>                                                      \
tmp<Field<Type> > wedgeFvPatchField<Type>::snGradTransformDiag()\
const                                                           \
{                                                               \
    return tmp<Field<Type> >                                    \
    (                                                           \
        new Field<Type>(this->size(), pTraits<Type>::zero)      \
    );                                                          \
}


#define doMakePatchTypeField(type, Type, args...)                           \
    makeVectorTensorNWedgeFuncDefs(type)                                    \
                                                                            \
    makePatchTypeField(fvPatch##Type##Field, wedgeFvPatch##Type##Field);


forAllVectorNTypes(doMakePatchTypeField)

forAllTensorNTypes(doMakePatchTypeField)

forAllDiagTensorNTypes(doMakePatchTypeField)

forAllSphericalTensorNTypes(doMakePatchTypeField)


#undef doMakePatchTypeField

#undef makeVectorTensorNWedgeFuncDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
