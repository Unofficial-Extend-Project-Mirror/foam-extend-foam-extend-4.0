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
    Foam::cyclicFaePatchField

Description
    Foam::cyclicFaePatchField

SourceFiles
    cyclicFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFaePatchField_H
#define cyclicFaePatchField_H

#include "coupledFaePatchField.H"
#include "cyclicFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFaePatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicFaePatchField
:
    public coupledFaePatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicFaPatch& cyclicPatch_;


public:

    //- Runtime type information
    TypeName(cyclicFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicFaePatchField onto a new patch
        cyclicFaePatchField
        (
            const cyclicFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        cyclicFaePatchField
        (
            const cyclicFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type> > clone() const
        {
            return tmp<faePatchField<Type> >
            (
                new cyclicFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicFaePatchField
        (
            const cyclicFaePatchField<Type>&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faePatchField<Type> > clone
        (
            const DimensionedField<Type, edgeMesh>& iF
        ) const
        {
            return tmp<faePatchField<Type> >
            (
                new cyclicFaePatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cyclicFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
