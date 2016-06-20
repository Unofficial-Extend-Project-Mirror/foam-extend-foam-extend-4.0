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

Class
    SymmetryPointPatchField

Description

SourceFiles
    SymmetryPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef SymmetryPointPatchField_H
#define SymmetryPointPatchField_H

#include "BasicSymmetryPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class SymmetryPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SymmetryPointPatch,
    template<class> class MatrixType,
    class Type
>
class SymmetryPointPatchField
:
    public BasicSymmetryPointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>
{

public:

    //- Runtime type information
    TypeName(SymmetryPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        SymmetryPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        SymmetryPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        SymmetryPointPatchField
        (
            const SymmetryPointPatchField
            <
                PatchField,
                Mesh,
                PointPatch,
                SymmetryPointPatch,
                MatrixType,
                Type
            >&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        SymmetryPointPatchField
        (
            const SymmetryPointPatchField
            <
                PatchField,
                Mesh,
                PointPatch,
                SymmetryPointPatch,
                MatrixType,
                Type
            >&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new SymmetryPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    SymmetryPointPatch,
                    MatrixType,
                    Type
                >(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SymmetryPointPatchField
        (
            const SymmetryPointPatchField
            <
                PatchField,
                Mesh,
                PointPatch,
                SymmetryPointPatch,
                MatrixType,
                Type
            >&,
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
                new SymmetryPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    SymmetryPointPatch,
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
#    include "SymmetryPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
