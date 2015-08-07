/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    elementPatchField

Description
    Abstract base class for FEM mesh element patch fields.

SourceFiles
    elementPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef elementPatchField_H
#define elementPatchField_H

#include "PointPatchFieldTemplate.H"
#include "tetPolyPatch.H"
#include "DummyMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
template<class> class DummyMatrix;
class elementMesh;

/*---------------------------------------------------------------------------*\
                       Class elementPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class elementPatchField
:
    public PointPatchField
    <
        elementPatchField,
        elementMesh,
        tetPolyPatch,
        DummyMatrix,
        Type
    >
{

public:

    typedef elementMesh Mesh;


    // Constructors

        //- Construct from patch and internal field
        elementPatchField
        (
            const tetPolyPatch&,
            const DimensionedField<Type, elementMesh>&
        );

        //- Construct as copy
        elementPatchField(const elementPatchField<Type>&);

        //- Construct and return a clone
        virtual autoPtr<elementPatchField<Type> > clone() const
        {
            return autoPtr<elementPatchField<Type> >
            (
                new elementPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        elementPatchField
        (
            const elementPatchField<Type>&,
            const DimensionedField<Type, elementMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<elementPatchField<Type> > clone
        (
            const DimensionedField<Type, elementMesh>& iF
        ) const
        {
            return autoPtr<elementPatchField<Type> >
            (
                new elementPatchField<Type>(*this, iF)
            );
        }


    // Destructor

        virtual ~elementPatchField()
        {}


    // Member operators

        virtual void operator=(const Type& t)
        {
            PointPatchField
            <
                Foam::elementPatchField,
                elementMesh,
                tetPolyPatch,
                Foam::DummyMatrix,
                Type
            >::operator=(t);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "elementPatchField.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
