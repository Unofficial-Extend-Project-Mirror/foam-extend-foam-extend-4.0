/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    cohesiveZoneFvPatchVectorField

Description

SourceFiles
    cohesiveZoneFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveZoneFvPatchVectorField_H
#define cohesiveZoneFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "cohesiveLaw.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cohesiveZoneFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Name of the displacement field
        const word UName_;

        //- Name of rheology model
        const word rheologyName_;

        //- Cohesive law
        cohesiveLaw* cohesiveLawPtr_;

        //- Under-relaxation factor
        scalar relaxationFactor_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("cohesiveZone");


    // Constructors

        //- Construct from patch and internal field
        cohesiveZoneFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cohesiveZoneFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cohesiveZoneFvPatchVectorField onto
        //  a new patch
        cohesiveZoneFvPatchVectorField
        (
            const cohesiveZoneFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new cohesiveZoneFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cohesiveZoneFvPatchVectorField
        (
            const cohesiveZoneFvPatchVectorField&,
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
                new cohesiveZoneFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

//         virtual ~cohesiveZoneFvPatchVectorField();


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return cohesive law
            const cohesiveLaw& law() const
            {
                return *cohesiveLawPtr_;
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
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
    

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
