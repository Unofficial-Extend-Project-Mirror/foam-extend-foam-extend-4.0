/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    planeContactDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation (stressedFoam).

SourceFiles
    planeContactDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef planeContactDisplacementFvPatchVectorField_H
#define planeContactDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class planeContactDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class planeContactDisplacementFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{

    // Private Data

    //- Plane initial position
    vector p0_;

    //- Plane normal
    vector n0_;

    //- Plane velocity in normal direction
    scalar U0_;

public:

    //- Runtime type information
    TypeName("planeContactDisplacement");


    // Constructors

        //- Construct from patch and internal field
        planeContactDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        planeContactDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  planeContactDisplacementFvPatchVectorField
        //  onto a new patch
        planeContactDisplacementFvPatchVectorField
        (
            const planeContactDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new planeContactDisplacementFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        planeContactDisplacementFvPatchVectorField
        (
            const planeContactDisplacementFvPatchVectorField&,
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
                new planeContactDisplacementFvPatchVectorField
                (
                    *this, iF
                )
            );
        }


    // Member functions

        // Access


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

            //- evaluate the patch
            virtual void evaluate(const Pstream::commsTypes);

            //- snGrad with non-orthogonal correction
            tmp<Foam::Field<vector> > snGrad() const;

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
