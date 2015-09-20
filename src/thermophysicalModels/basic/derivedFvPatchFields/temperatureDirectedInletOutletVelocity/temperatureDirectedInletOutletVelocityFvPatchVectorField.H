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
    Foam::temperatureDirectedInletOutletVelocityFvPatchVectorField

Description
    Foam::temperatureDirectedInletOutletVelocityFvPatchVectorField

SourceFiles
    temperatureDirectedInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDirectedInletOutletVelocityFvPatchVectorField_H
#define temperatureDirectedInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class temperatureDirectedInletOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class temperatureDirectedInletOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Name of flux field
        word phiName_;

        //- Name of static Temperature field
        word TName_;

        //- Total temperature field
        scalarField T0_;

        //- Inlet direction
        vectorField inletDir_;

        //- Is the supplied inlet value in cartesian or
        //  cylindrical coordinates?
        Switch cylindricalCCS_;

        //- Angular velocity of the frame
        vector omega_;


public:

    //- Runtime type information
    TypeName("temperatureDirectedInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        temperatureDirectedInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureDirectedInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  temperatureDirectedInletOutletVelocityFvPatchVectorField
        //  onto a new patch
        temperatureDirectedInletOutletVelocityFvPatchVectorField
        (
            const temperatureDirectedInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperatureDirectedInletOutletVelocityFvPatchVectorField
        (
            const temperatureDirectedInletOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new temperatureDirectedInletOutletVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        temperatureDirectedInletOutletVelocityFvPatchVectorField
        (
            const temperatureDirectedInletOutletVelocityFvPatchVectorField&,
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
                new temperatureDirectedInletOutletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


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
