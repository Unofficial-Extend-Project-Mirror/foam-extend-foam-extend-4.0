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
    Foam::edgeNormalFixedValueFaPatchVectorField

Description
    Edge normal fixed value vector field finite area boundary condition

    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    edgeNormalFixedValueFaPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef edgeNormalFixedValueFaPatchVectorField_H
#define edgeNormalFixedValueFaPatchVectorField_H

#include "faPatchFields.H"
#include "fixedValueFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class edgeNormalFixedValueFaPatch Declaration
\*---------------------------------------------------------------------------*/

class edgeNormalFixedValueFaPatchVectorField
:
    public fixedValueFaPatchVectorField
{
    // Private data

        //- Surface-normal velocity value
        scalarField refValue_;


public:

    //- Runtime type information
    TypeName("edgeNormalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  edgeNormalFixedValueFaPatchVectorField
        //  onto a new patch
        edgeNormalFixedValueFaPatchVectorField
        (
            const edgeNormalFixedValueFaPatchVectorField&,
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        edgeNormalFixedValueFaPatchVectorField
        (
            const edgeNormalFixedValueFaPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<faPatchVectorField> clone() const
        {
            return tmp<faPatchVectorField>
            (
                new edgeNormalFixedValueFaPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        edgeNormalFixedValueFaPatchVectorField
        (
            const edgeNormalFixedValueFaPatchVectorField&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchVectorField> clone
        (
            const DimensionedField<vector, areaMesh>& iF
        ) const
        {
            return tmp<faPatchVectorField>
            (
                new edgeNormalFixedValueFaPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchVectorField&,
                const labelList&
            );


    // Member functions

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
