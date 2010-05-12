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
    emptyFaPatchField

Description

SourceFiles
    emptyFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef emptyFaPatchField_H
#define emptyFaPatchField_H

#include "faPatchField.H"
#include "emptyFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class emptyFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class emptyFaPatchField
:
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(emptyFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        emptyFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        emptyFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given emptyFaPatchField onto a new patch
        emptyFaPatchField
        (
            const emptyFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        emptyFaPatchField
        (
            const emptyFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >
            (
                new emptyFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        emptyFaPatchField
        (
            const emptyFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type> >
            (
                new emptyFaPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            )
            {}

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
                const labelList&
            )
            {}


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  This only checks to see the case is actually 1D or 2D
            //  for which this boundary condition is valid
            void updateCoeffs();


            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type> > gradientInternalCoeffs() const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type> > gradientBoundaryCoeffs() const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "emptyFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
