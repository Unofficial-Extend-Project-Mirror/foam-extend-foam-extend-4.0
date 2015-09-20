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
    Foam::incompressible::RASModels::RWallFunctionFvPatchSymmTensorField

Description
    Boundary condition for R when using wall functions
    - set and updated from the epsilon wall function

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    RWallFunctionFvPatchSymmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef RWallFunctionFvPatchSymmTensorField_H
#define RWallFunctionFvPatchSymmTensorField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
             Class RWallFunctionFvPatchSymmTensorField Declaration
\*---------------------------------------------------------------------------*/

class RWallFunctionFvPatchSymmTensorField
:
    public zeroGradientFvPatchSymmTensorField
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of turbulent viscosity field
        word nutName_;


    // Private member functions

        //- Check the type of the patch
        void checkType();

        //- Set value in internal field
        void setInInternalField(const symmTensorField& f) const;


public:

    //- Runtime type information
    TypeName("RWallFunction");


    // Constructors

        //- Construct from patch and internal field
        RWallFunctionFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RWallFunctionFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  RWallFunctionFvPatchSymmTensorField
        //  onto a new patch
        RWallFunctionFvPatchSymmTensorField
        (
            const RWallFunctionFvPatchSymmTensorField&,
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        RWallFunctionFvPatchSymmTensorField
        (
            const RWallFunctionFvPatchSymmTensorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchSymmTensorField> clone() const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new RWallFunctionFvPatchSymmTensorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        RWallFunctionFvPatchSymmTensorField
        (
            const RWallFunctionFvPatchSymmTensorField&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchSymmTensorField> clone
        (
            const DimensionedField<symmTensor, volMesh>& iF
        ) const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new RWallFunctionFvPatchSymmTensorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
