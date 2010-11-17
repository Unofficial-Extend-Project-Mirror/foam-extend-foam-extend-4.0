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
    Foam::symmetryFaePatchField

Description
    Foam::symmetryFaePatchField

SourceFiles
    symmetryFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryFaePatchField_H
#define symmetryFaePatchField_H

#include "faePatchField.H"
#include "symmetryFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class symmetryFaePatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class symmetryFaePatchField
:
    public faePatchField<Type>
{

public:

    //- Runtime type information
    TypeName(symmetryFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        symmetryFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetryFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given symmetryFaePatchField onto a new patch
        symmetryFaePatchField
        (
            const symmetryFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        symmetryFaePatchField
        (
            const symmetryFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type> > clone() const
        {
            return tmp<faePatchField<Type> >
            (
                new symmetryFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        symmetryFaePatchField
        (
            const symmetryFaePatchField<Type>&,
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
                new symmetryFaePatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "symmetryFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
