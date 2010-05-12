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
    transformFaPatchField

Description

SourceFiles
    transformFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef transformFaPatchField_H
#define transformFaPatchField_H

#include "faPatchField.H"
#include "areaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class transformFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class transformFaPatchField
:
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("transform");


    // Constructors

        //- Construct from patch and internal field
        transformFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transformFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given transformFaPatchField<Type>
        //  onto a new patch
        transformFaPatchField
        (
            const transformFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        transformFaPatchField
        (
            const transformFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>&
        ) const = 0;


    // Member functions

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type> > snGrad() const = 0;

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type> > snGradTransformDiag() const = 0;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


    // Member operators

        virtual void operator=(const faPatchField<Type>&);
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField > transformFaPatchField<scalar>::gradientInternalCoeffs() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "transformFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
