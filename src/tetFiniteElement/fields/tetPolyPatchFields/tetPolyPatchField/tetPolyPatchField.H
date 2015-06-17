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
    tetPolyPatchField

Description
    Abstract base class for tetrahedral mesh patch fields. Note the special
    mechanism at work here: the field itself holds no values, as the point
    values belong to the internal field. However, the field will create a
    list of values if required.

SourceFiles
    tetPolyPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchField_H
#define tetPolyPatchField_H

#include "PointPatchFieldTemplate.H"
#include "tetPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
template<class> class tetFemMatrix;
class tetPointMesh;

/*---------------------------------------------------------------------------*\
                       Class tetPolyPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class tetPolyPatchField
:
    public PointPatchField
    <
        tetPolyPatchField,
        tetPointMesh,
        tetPolyPatch,
        tetFemMatrix,
        Type
    >
{

public:

    // Constructors

        //- Construct from patch and internal field
        tetPolyPatchField
        (
            const tetPolyPatch&,
            const DimensionedField<Type, tetPointMesh>&
        );

        //- Construct as copy
        tetPolyPatchField(const tetPolyPatchField<Type>&);

        //- Construct and return a clone
        virtual autoPtr<tetPolyPatchField<Type> > clone() const
        {
            return autoPtr<tetPolyPatchField<Type> >
            (
                new tetPolyPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tetPolyPatchField
        (
            const tetPolyPatchField<Type>&,
            const DimensionedField<Type, tetPointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<tetPolyPatchField<Type> > clone
        (
            const DimensionedField<Type, tetPointMesh>& iF
        ) const
        {
            return autoPtr<tetPolyPatchField<Type> >
            (
                new tetPolyPatchField<Type>(*this, iF)
            );
        }


    // Destructor

        virtual ~tetPolyPatchField()
        {}


    // Member operators

        virtual void operator=(const Type& t)
        {
            PointPatchField
            <
                Foam::tetPolyPatchField,
                tetPointMesh,
                tetPolyPatch,
                Foam::tetFemMatrix,
                Type
            >::operator=(t);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetPolyPatchField.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

