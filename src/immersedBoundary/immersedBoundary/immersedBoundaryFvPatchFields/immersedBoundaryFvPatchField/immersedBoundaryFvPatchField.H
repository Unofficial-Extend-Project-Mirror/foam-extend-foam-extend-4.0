/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::immersedBoundaryFvPatchField

Description
    Constrained immersedBoundaryFvPatchField of calculated functionality
    with appropriate resizing on topo change or motion of the immersed surface.
    The patch field shall take the values of intersected cells.

    Note: This patch field is used instead of the calculated type in the
    immersedBoundaryFvPatch.  It does not support discretisation.
    For matrix operations, use the mixedIbFvPatchField.

Author
    Hrvoje Jasak

SourceFiles
    immersedBoundaryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvPatchField_H
#define immersedBoundaryFvPatchField_H

#include "fvPatchFields.H"
#include "immersedBoundaryFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryFvPatchField
:
    public fvPatchField<Type>,
    public immersedBoundaryFieldBase<Type>
{
public:

    //- Runtime type information: constrained type
    TypeName(immersedBoundaryFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given immersedBoundaryFvPatchField
        //  onto a new patch
        immersedBoundaryFvPatchField
        (
            const immersedBoundaryFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryFvPatchField
        (
            const immersedBoundaryFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryFvPatchField
        (
            const immersedBoundaryFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~immersedBoundaryFvPatchField()
    {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );

            //- Update on mesh motion
            virtual void updateOnMotion();


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
