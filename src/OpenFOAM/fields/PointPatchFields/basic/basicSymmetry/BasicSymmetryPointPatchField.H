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
    BasicSymmetryPointPatchField

Description

SourceFiles
    BasicSymmetryPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef BasicSymmetryPointPatchField_H
#define BasicSymmetryPointPatchField_H

#include "ValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class BasicSymmetryPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
class BasicSymmetryPointPatchField
:
    public ValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
{
    // Private memeber functions

        //- Calculate implicitness matrix
        Type fixFlagMatrix(const vector&) const;


public:

    // Constructors

        //- Construct from patch and internal field
        BasicSymmetryPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        BasicSymmetryPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        BasicSymmetryPointPatchField
        (
            const BasicSymmetryPointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        BasicSymmetryPointPatchField
        (
            const BasicSymmetryPointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new BasicSymmetryPointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        BasicSymmetryPointPatchField
        (
            const BasicSymmetryPointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&,
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
                new BasicSymmetryPointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>(*this, iF)
            );
        }


    // Member functions

        //- Update the patch field
        virtual void updateBoundaryField();

        //- Set boundary condition to matrix
        virtual void setBoundaryCondition
        (
            Map<typename MatrixType<Type>::ConstraintType>&
        ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "BasicSymmetryPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
