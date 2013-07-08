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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidTractionFreeFvPatchVectorField

Description
    Traction free boundary condition for solid foam solvers.
    Suitable for for total strain or incremental strain total or
    updated Lagrangian approaches.
    Contains thermal-elastic, plastic, non-linear terms if required.
    And non-orthogonal correction.

Author
   Philip Cardiff

SourceFiles
    solidTractionFreeFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionFreeFvPatchVectorField_H
#define solidTractionFreeFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidTractionFreeFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionFreeFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- non-linear solver options
        enum nonLinearType
	  {
	    OFF,
	    UPDATED_LAGRANGIAN,
	    TOTAL_LAGRANGIAN
	  };
  
        static const NamedEnum<nonLinearType, 3> nonLinearNames_;

        //- if it is a non linear solver
        nonLinearType nonLinear_;

public:

    //- Runtime type information
    TypeName("solidTractionFree");


    // Constructors

        //- Construct from patch and internal field
        solidTractionFreeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionFreeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidTractionFreeFvPatchVectorField onto a new patch
        solidTractionFreeFvPatchVectorField
        (
            const solidTractionFreeFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidTractionFreeFvPatchVectorField
        (
            const solidTractionFreeFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidTractionFreeFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidTractionFreeFvPatchVectorField
        (
            const solidTractionFreeFvPatchVectorField&,
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
                new solidTractionFreeFvPatchVectorField(*this, iF)
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


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
