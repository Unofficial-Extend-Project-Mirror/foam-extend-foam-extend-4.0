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
    EmptyPointPatchField

Description

SourceFiles
    EmptyPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef EmptyPointPatchField_H
#define EmptyPointPatchField_H

#include "PointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class EmptyPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class EmptyPointPatch,
    template<class> class MatrixType,
    class Type
>
class EmptyPointPatchField
:
    public PatchField<Type>
{

public:

    //- Runtime type information
    TypeName(EmptyPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        EmptyPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        EmptyPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        EmptyPointPatchField
        (
            const EmptyPointPatchField
            <PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        EmptyPointPatchField
        (
            const EmptyPointPatchField
            <PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new EmptyPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    EmptyPointPatch,
                    MatrixType,
                    Type
                >(*this)
            );
        }

        //- Construct as copy setting internal field reference
        EmptyPointPatchField
        (
            const EmptyPointPatchField
            <PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const DimensionedField<Type, Mesh>& iF
        ) const
        {
            return autoPtr<PatchField<Type> >
            (
                new EmptyPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    EmptyPointPatch,
                    MatrixType,
                    Type
                >(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "EmptyPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
