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
    pRveTractionDisplacementIncrementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation (stressFoam).

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    pRveTractionDisplacementIncrementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pRveTractionDisplacementIncrementFvPatchVectorField_H
#define pRveTractionDisplacementIncrementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class pRveTractionDisplacementIncrementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pRveTractionDisplacementIncrementFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;


public:

    //- Runtime type information
    TypeName("pRveTractionDisplacementIncrement");


    // Constructors

        //- Construct from patch and internal field
        pRveTractionDisplacementIncrementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pRveTractionDisplacementIncrementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pRveTractionDisplacementIncrementFvPatchVectorField
        //  onto a new patch
        pRveTractionDisplacementIncrementFvPatchVectorField
        (
            const pRveTractionDisplacementIncrementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pRveTractionDisplacementIncrementFvPatchVectorField
        (
            const pRveTractionDisplacementIncrementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pRveTractionDisplacementIncrementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pRveTractionDisplacementIncrementFvPatchVectorField
        (
            const pRveTractionDisplacementIncrementFvPatchVectorField&,
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
                new pRveTractionDisplacementIncrementFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual  scalarField& pressure()
            {
                return pressure_;
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
