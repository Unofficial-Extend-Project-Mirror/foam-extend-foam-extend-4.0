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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF),   // Do not read mixed data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "immersedBoundaryFvPatchField<Type>::"
            "immersedBoundaryFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(p, iF),  // Do not map base data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatchField<Type>::"
            "immersedBoundaryFvPatchField\n"
            "(\n"
            "    const immersedBoundaryFvPatchField<Type>&,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    // On creation of the mapped field, the internal field is dummy and
    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    Field<Type>::operator=(Field<Type>(p.size(), pTraits<Type>::zero));
}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    ibPatch_(ptf.ibPatch())
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    ibPatch_(ptf.ibPatch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void immersedBoundaryFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Use internal values
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{
    // Base fields do not rmap
    this->setSize(this->patch().size(), pTraits<Type>::zero);
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    // Use internal values
    Field<Type>::operator=(this->patchInternalField());

    fvPatchField<Type>::evaluate();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::write(Ostream& os) const
{
    // Resolve post-processing issues.  HJ, 1/Dec/2017
    fvPatchField<Type>::write(os);

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    // Write VTK on master only
    if (Pstream::master())
    {
        // Add parallel reduction of all faces and data to proc 0
        // and write the whola patch together

        // Write immersed boundary data as a vtk file
        autoPtr<surfaceWriter> writerPtr = surfaceWriter::New("vtk");

        // Get the intersected patch
        const standAlonePatch& ts = ibPatch_.ibPolyPatch().ibPatch();

        writerPtr->write
        (
            this->dimensionedInternalField().path(),
            ibPatch_.name(),
            ts.points(),
            ts,
            this->dimensionedInternalField().name(),
            *this,
            false, // FACE_DATA
            false  // verbose
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
