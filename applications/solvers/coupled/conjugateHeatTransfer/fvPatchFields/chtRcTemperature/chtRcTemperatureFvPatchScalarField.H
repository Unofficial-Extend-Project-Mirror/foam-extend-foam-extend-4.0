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
    chtRcTemperatureFvPatchScalarField

Description
    Region couple patch field taking radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    chtRcTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcTemperatureFvPatchScalarField_H
#define chtRcTemperatureFvPatchScalarField_H

#include "regionCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class chtRcTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcTemperatureFvPatchScalarField
:
    public regionCouplingFvPatchScalarField
{
    // Private data

        //- Name of thermal diffusivity field
        word kName_;

        //- Take radiation into account for this patch
        bool radiation_;


public:

    //- Runtime type information
    TypeName("chtRcTemperature");


    // Constructors

        //- Construct from patch and internal field
        chtRcTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chtRcTemperatureFvPatchScalarField onto a new patch
        chtRcTemperatureFvPatchScalarField
        (
            const chtRcTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcTemperatureFvPatchScalarField
        (
            const chtRcTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            //- Return radiation
            bool radiation() const
            {
                return radiation_;
            }

            //- Return shadow patch field
            virtual const chtRcTemperatureFvPatchScalarField&
                shadowPatchField() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Calculate source term
            virtual tmp<scalarField> source() const;

            //- Manipulate matrix
            virtual void manipulateMatrix(fvScalarMatrix& matrix);

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
