/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::decayingTurbulenceFvPatchVectorField

Description
    Decaying turbulence generator for LES, using mapped data for profiles

SourceFiles
    decayingTurbulenceFvPatchVectorField.C

Author
    Prof. Nikolai Kornev, Uni Rostock
    Clean-up by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef decayingTurbulenceFvPatchVectorField_H
#define decayingTurbulenceFvPatchVectorField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "SLList.H"

#include "decayingVorton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class decayingTurbulenceFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class decayingTurbulenceFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Length-scale field
        scalarField LField_;

        //- Reference velocity field
        vectorField refField_;

        //- Reference Reynolds stress field
        symmTensorField RField_;

        //- Field
        tensorField Lund_;

        //- Reynolds stress field
        symmTensorField R_;

        //- Direction
        label direction_;

        //- Vortex overlap (defaults to 0.3)
        scalar vortexOverlap_;

        //- Number of vortexes (defaults to 10)
        label nVortexes_;

        //- Time index for update
        label curTimeIndex_;

        //- Random number generator
        Random rndGen_;

        //- Decaying vortons.  In parallel, stored only on master processor
        SLList<decayingVorton> vortons_;

        //- Index
        label index_;


    // Private member functions

        //- Update vortons
        void updateVortons();

        //- Lspot
        static inline scalar lSpot(const scalar l)
        {
            return 3*l;
        }


public:

    //- Runtime type information
    TypeName("decayingTurbulenceInflowGenerator");


    // Constructors

        //- Construct from patch and internal field
        decayingTurbulenceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        decayingTurbulenceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        decayingTurbulenceFvPatchVectorField
        (
            const decayingTurbulenceFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        decayingTurbulenceFvPatchVectorField
        (
            const decayingTurbulenceFvPatchVectorField&
        );

        //- Construct and return a clone setting internal field reference
        decayingTurbulenceFvPatchVectorField
        (
            const decayingTurbulenceFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new decayingTurbulenceFvPatchVectorField(*this)
            );
        }


        //- Construct as copy setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new decayingTurbulenceFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        // Evaluation functions

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
