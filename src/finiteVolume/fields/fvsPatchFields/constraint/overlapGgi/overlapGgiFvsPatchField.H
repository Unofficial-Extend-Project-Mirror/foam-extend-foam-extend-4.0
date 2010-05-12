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

Class
    overlapGgiFvsPatchField

Description
    Partial overlap generalized grid interface patch field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    overlapGgiFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiFvsPatchField_H
#define overlapGgiFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "overlapGgiFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class overlapGgiFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class overlapGgiFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the overlapGgi patch
        const overlapGgiFvPatch& overlapGgiPatch_;


public:

    //- Runtime type information
    TypeName(overlapGgiFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        overlapGgiFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        overlapGgiFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given overlapGgiFvsPatchField onto a new patch
        overlapGgiFvsPatchField
        (
            const overlapGgiFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        overlapGgiFvsPatchField
        (
            const overlapGgiFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new overlapGgiFvsPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new overlapGgiFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "overlapGgiFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
