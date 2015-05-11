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
    Foam::immersedBoundaryFvsPatchField

Description
    Foam::immersedBoundaryFvsPatchField

Author
    Hrvoje Jasak

SourceFiles
    immersedBoundaryFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvsPatchField_H
#define immersedBoundaryFvsPatchField_H

#include "fvsPatchField.H"
#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class immersedBoundaryFvPatchField;

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryFvsPatchField
:
    public fvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const immersedBoundaryFvPatch& ibPatch_;


public:

    //- Runtime type information
    TypeName(immersedBoundaryFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given immersedBoundaryFvsPatchField
        //  onto a new patch
        immersedBoundaryFvsPatchField
        (
            const immersedBoundaryFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryFvsPatchField
        (
            const immersedBoundaryFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new immersedBoundaryFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryFvsPatchField
        (
            const immersedBoundaryFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new immersedBoundaryFvsPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~immersedBoundaryFvsPatchField()
    {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            )
            {}

            //- Reverse map the given fvsPatchField onto this fvsPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            )
            {}


    // Member functions

        //- Return reference to immersed boundary patch
        const immersedBoundaryFvPatch& ibPatch() const
        {
            return ibPatch_;
        }

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

//         virtual void operator=(const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
