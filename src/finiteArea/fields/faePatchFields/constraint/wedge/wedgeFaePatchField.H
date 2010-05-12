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
    wedgeFaePatchField

Description

SourceFiles
    wedgeFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFaePatchField_H
#define wedgeFaePatchField_H

#include "faePatchField.H"
#include "wedgeFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class wedgeFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class wedgeFaePatchField
:
    public faePatchField<Type>
{

public:

    //- Runtime type information
    TypeName(wedgeFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        wedgeFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wedgeFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wedgeFaePatchField onto a new patch
        wedgeFaePatchField
        (
            const wedgeFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        wedgeFaePatchField
        (
            const wedgeFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type> > clone() const
        {
            return tmp<faePatchField<Type> >
            (
                new wedgeFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wedgeFaePatchField
        (
            const wedgeFaePatchField<Type>&,
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
                new wedgeFaePatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "wedgeFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
