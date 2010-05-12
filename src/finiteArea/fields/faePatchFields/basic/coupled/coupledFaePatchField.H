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
    coupledFaePatchField

Description

SourceFiles
    coupledFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFaePatchField_H
#define coupledFaePatchField_H

#include "faePatchField.H"
#include "coupledFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coupledFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledFaePatchField
:
    public faePatchField<Type>
{

public:

    //- Runtime type information
    TypeName(coupledFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        coupledFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch and internal field and patch field
        coupledFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        coupledFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given coupledFaePatchField onto a new patch
        coupledFaePatchField
        (
            const coupledFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        coupledFaePatchField
        (
            const coupledFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        coupledFaePatchField
        (
            const coupledFaePatchField<Type>&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type> > clone
        (
            const DimensionedField<Type, edgeMesh>&
        ) const = 0;


    // Member functions

        // Access

            //- Return true if this patch field is derived from
            //  coupledFaePatchField<Type>.
            virtual bool coupled() const
            {
                return true;
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "coupledFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
