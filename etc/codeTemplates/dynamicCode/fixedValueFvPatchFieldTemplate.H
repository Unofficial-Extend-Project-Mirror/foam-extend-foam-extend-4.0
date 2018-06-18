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

Description
    Template for use with dynamic code generation of a
    fixedValue fvPatchField.

    - without state

SourceFiles
    fixedValueFvPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueFvPatchTemplate${FieldType}_H
#define fixedValueFvPatchTemplate${FieldType}_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated FixedValueFvPatch
\*---------------------------------------------------------------------------*/

class ${typeName}FixedValueFvPatch${FieldType}
:
    public fixedValueFvPatchField<${TemplateType}>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const ${typeName}FixedValueFvPatch${FieldType}&,
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const ${typeName}FixedValueFvPatch${FieldType}&
        );

        //- Construct and return a clone
        virtual tmp< fvPatch${FieldType} > clone() const
        {
            return tmp< fvPatch${FieldType} >
            (
                new ${typeName}FixedValueFvPatch${FieldType}(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const ${typeName}FixedValueFvPatch${FieldType}&,
            const DimensionedField<${TemplateType}, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp< fvPatch${FieldType} > clone
        (
            const DimensionedField<${TemplateType}, volMesh>& iF
        ) const
        {
            return tmp< fvPatch${FieldType} >
            (
                new ${typeName}FixedValueFvPatch${FieldType}(*this, iF)
            );
        }


    //- Destructor
    virtual ~${typeName}FixedValueFvPatch${FieldType}();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
