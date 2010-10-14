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

#ifndef wedgeFvPatchVectorNFields_H
#define wedgeFvPatchVectorNFields_H

#include "volVectorNFields.H"
#include "fvPatchVectorNFields.H"

#include "wedgeFvPatchVectorNFieldsFwd.H"
#include "wedgeFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define doMakeVectorTensorNWedgeFuncDefs(type, Type, args...)       \
template<>                                                          \
tmp<Field<type > > wedgeFvPatchField<type >::snGrad() const;        \
                                                                    \
template<>                                                          \
void wedgeFvPatchField<type >::evaluate(                            \
    const Pstream::commsTypes commsType                             \
);                                                                  \
                                                                    \
template<>                                                          \
tmp<Field<type > > wedgeFvPatchField<type >::snGradTransformDiag()  \
const;

forAllVectorNTypes(doMakeVectorTensorNWedgeFuncDefs)

forAllTensorNTypes(doMakeVectorTensorNWedgeFuncDefs)

forAllDiagTensorNTypes(doMakeVectorTensorNWedgeFuncDefs)

forAllSphericalTensorNTypes(doMakeVectorTensorNWedgeFuncDefs)

#undef doMakeVectorTensorNWedgeFuncDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
