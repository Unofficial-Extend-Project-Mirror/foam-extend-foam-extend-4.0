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
    WedgePointPatchField

Description
    Wedge front and back plane patch field

SourceFiles
    WedgePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef WedgePointPatchField_H
#define WedgePointPatchField_H

#include "ZeroGradientPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class WedgePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
class WedgePointPatchField
:
    public PatchField<Type>
{

public:

    //- Runtime type information
    TypeName(WedgePointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        WedgePointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        WedgePointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        WedgePointPatchField
        (
            const WedgePointPatchField
            <PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        WedgePointPatchField
        (
            const WedgePointPatchField
            <PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new WedgePointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    WedgePointPatch,
                    MatrixType,
                    Type
                >(*this)
            );
        }

        //- Construct as copy setting internal field reference
        WedgePointPatchField
        (
            const WedgePointPatchField
            <PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>&,
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
                new WedgePointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    WedgePointPatch,
                    MatrixType,
                    Type
                >(*this, iF)
            );
        }


    // Member functions

        //- Update the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::blocking
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "WedgePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
