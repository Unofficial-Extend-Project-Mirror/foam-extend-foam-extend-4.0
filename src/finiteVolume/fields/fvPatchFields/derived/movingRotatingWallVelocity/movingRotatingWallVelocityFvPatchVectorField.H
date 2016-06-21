/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    movingRotatingWallVelocityFvPatchVectorField

Description
    movingRotatingWallVelocityFvPatchVectorField

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    movingRotatingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingRotatingWallVelocityFvPatchVectorField_H
#define movingRotatingWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class movingRotatingWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingRotatingWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Centre of rotation
        point centre_;

        //- Rotation axis
        vector axis_;

        //- Rotational speed in rpm
        scalar rpm_;


public:

    //- Runtime type information
    TypeName("movingRotatingWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        movingRotatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from components
        movingRotatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const point& centre,
            const vector& axis,
            const scalar rpm
        );

        //- Construct from patch, internal field and dictionary
        movingRotatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingRotatingWallVelocityFvPatchVectorField
        //  onto a new patch
        movingRotatingWallVelocityFvPatchVectorField
        (
            const movingRotatingWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingRotatingWallVelocityFvPatchVectorField
        (
            const movingRotatingWallVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingRotatingWallVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingRotatingWallVelocityFvPatchVectorField
        (
            const movingRotatingWallVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new movingRotatingWallVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access

            //- Return centre of rotation
            const point& centre() const
            {
                return centre_;
            }

            //- Return axis of rotation
            const vector& axis() const
            {
                return axis_;
            }

            //- Return rpm
            scalar rpm() const
            {
                return rpm_;
            }


        // Edit

            //- Return access to centre of rotation
            point& centre()
            {
                return centre_;
            }

            //- Return access to axis of rotation
            vector& axis()
            {
                return axis_;
            }

            //- Return access to rpm
            scalar& rpm()
            {
                return rpm_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
