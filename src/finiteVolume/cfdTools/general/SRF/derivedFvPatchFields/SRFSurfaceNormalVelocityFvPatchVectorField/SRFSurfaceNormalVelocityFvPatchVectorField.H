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
    Foam::SRFSurfaceNormalVelocityFvPatchVectorField

Description
    Velocity patch, specifying surface-normal velocity to be used with the
    SRF model.  Surface normal vector boundary condition is set by its
    magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    SRFSurfaceNormalVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef SRFSurfaceNormalVelocityFvPatchVectorField_H
#define SRFSurfaceNormalVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class SRFSurfaceNormalVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class SRFSurfaceNormalVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Is the supplied inlet value relative to the SRF
        Switch relative_;

        //- Surface-normal velocity value
        scalarField refValue_;


public:

    //- Runtime type information
    TypeName("SRFSurfaceNormalVelocity");


    // Constructors

        //- Construct from patch and internal field
        SRFSurfaceNormalVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SRFSurfaceNormalVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  SRFSurfaceNormalVelocityFvPatchVectorField onto a new patch
        SRFSurfaceNormalVelocityFvPatchVectorField
        (
            const SRFSurfaceNormalVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        SRFSurfaceNormalVelocityFvPatchVectorField
        (
            const SRFSurfaceNormalVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new SRFSurfaceNormalVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SRFSurfaceNormalVelocityFvPatchVectorField
        (
            const SRFSurfaceNormalVelocityFvPatchVectorField&,
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
                new SRFSurfaceNormalVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return const access to the relative flag
            const Switch& relative() const
            {
                return relative_;
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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
