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
    pointPatchField

Description
    Abstract base class for point patch fields. Note the special
    mechanism at work here: the field itself holds no values, as the point
    values belong to the internal field. However, the field will create a
    list of values if required.

SourceFiles
    pointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatchField_H
#define pointPatchField_H

#include "PointPatchField.H"
#include "pointPatch.H"
#include "DummyMatrix.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class pointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointPatchField
:
    public PointPatchField
    <
        pointPatchField,
        pointMesh,
        pointPatch,
        DummyMatrix,
        Type
    >
{

public:

    // Constructors

        //- Construct from patch and internal field
        pointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct as copy
        pointPatchField(const pointPatchField<Type>&);

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type> > clone() const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new pointPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pointPatchField
        (
            const pointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type> > clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type> >
            (
                new pointPatchField<Type>(*this, iF)
            );
        }


    // Destructor

        virtual ~pointPatchField()
        {}


    // Member operators

        virtual void operator=(const Type& t)
        {
            PointPatchField
            <
                Foam::pointPatchField,
                pointMesh,
                pointPatch,
                DummyMatrix,
                Type
            >::operator=(t);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointPatchField.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

