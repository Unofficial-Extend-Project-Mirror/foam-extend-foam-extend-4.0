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
    Foam::solidSymmetryFvPatchVectorField

Description
    symmetry plane with non-orthogonal correction
    same as symmetryDisplacement

SourceFiles
    solidSymmetryFvPatchVectorField.C

Author
    Zeljko Tukovic FSB Zagreb

\*---------------------------------------------------------------------------*/

#ifndef solidSymmetryFvPatchVectorField_H
#define solidSymmetryFvPatchVectorField_H

#include "symmetryFvPatchFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solidSymmetryFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidSymmetryFvPatchVectorField
:
    public symmetryFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("solidSymmetry");


    // Constructors

        //- Construct from patch and internal field
        solidSymmetryFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidSymmetryFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidSymmetryFvPatchVectorField
        //  onto a new patch
        solidSymmetryFvPatchVectorField
        (
            const solidSymmetryFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidSymmetryFvPatchVectorField
        (
            const solidSymmetryFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidSymmetryFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidSymmetryFvPatchVectorField
        (
            const solidSymmetryFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidSymmetryFvPatchVectorField(*this, iF)
            );
        }

        tmp<Field<vector> > snGrad() const;

        void evaluate(const Pstream::commsTypes);

        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
