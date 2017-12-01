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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatchField.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void immersedBoundaryFvPatchField<Type>::updateIbValues()
{
    // Interpolate the values form tri surface using nearest triangle
    const labelList& nt = ibPatch_.ibPolyPatch().nearestTri();

    this->refValue() = Field<Type>(triValue_, nt);
    this->refGrad() = Field<Type>(triGrad_, nt);
    this->valueFraction() = scalarField(triValueFraction_, nt);
    Info<< "this->refValue(): " << this->refValue() << endl;
    mixedFvPatchField<Type>::evaluate();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::setDeadValues()
{
    // Fix the value in dead cells
    if (setDeadValue_)
    {
        const labelList& dc = ibPatch_.ibPolyPatch().deadCells();

        // Get non-const access to internal field
        Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

        forAll (dc, dcI)
        {
            psiI[dc[dcI]] = deadValue_;
        }
    }
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::correctDiag
(
    fvMatrix<Type>& eqn
) const
{
    scalarField& Diag = eqn.diag();

    const labelList& deadCells = ibPatch_.ibPolyPatch().deadCells();

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    if (deadCells.size() < Diag.size())
    {
        liveDiag = gSumMag(Diag)/(Diag.size() - deadCells.size());

        // Correct for sign
        liveDiag *= sign(gMax(Diag));
    }

    forAll (deadCells, cellI)
    {
        if (mag(Diag[deadCells[cellI]]) < SMALL)
        {
            Diag[deadCells[cellI]] = liveDiag;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    triValue_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    triGrad_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    triValueFraction_(false),
    setDeadValue_(false),
    deadValue_(pTraits<Type>::zero)
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),   // Do not read mixed data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    triValue_("triValue", dict, ibPatch_.ibMesh().size()),
    triGrad_("triGradient", dict, ibPatch_.ibMesh().size()),
    triValueFraction_("triValueFraction", dict, ibPatch_.ibMesh().size()),
    setDeadValue_(dict.lookup("setDeadValue")),
    deadValue_(pTraits<Type>(dict.lookup("deadValue")))
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

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();
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
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    triValue_(ptf.triValue()),
    triGrad_(ptf.triGrad()),
    triValueFraction_(ptf.triValueFraction()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
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

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();
}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    ibPatch_(ptf.ibPatch()),
    triValue_(ptf.triValue()),
    triGrad_(ptf.triGrad()),
    triValueFraction_(ptf.triValueFraction()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    ibPatch_(ptf.ibPatch()),
    triValue_(ptf.triValue()),
    triGrad_(ptf.triGrad()),
    triValueFraction_(ptf.triValueFraction()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void immersedBoundaryFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Base fields do not map: re-interpolate them from tri data
    this->updateIbValues();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{
    // Base fields do not rmap: re-interpolate them from tri data

    const immersedBoundaryFvPatchField<Type>& mptf =
        refCast<const immersedBoundaryFvPatchField<Type> >(ptf);

    // Set rmap tri data
    triValue_ = mptf.triValue_;
    triGrad_ = mptf.triGrad_;
    triValueFraction_ = mptf.triValueFraction_;

    this->updateIbValues();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    // Set dead value
    this->setDeadValues();

    // Evaluate mixed condition
    mixedFvPatchField<Type>::evaluate();
}


// template<class Type>
// void immersedBoundaryFvPatchField<Type>::manipulateMatrix
// (
//     fvMatrix<Type>& eqn
// )
// {
//     // Build matrix diagonal for cells where it is missing
//     this->correctDiag(eqn);

//     // Set values in IB cells
//     Field<Type> polyPsi(eqn.psi(), ibPatch_.ibCells());
//     eqn.setValues(ibPatch_.ibCells(), polyPsi);

//     // Correct equation for dead cells
//     Field<Type> deadCellsPsi
//     (
//         ibPatch_.deadCells().size(),
//         deadValue_
//     );
//     eqn.setValues(ibPatch_.deadCells(), deadCellsPsi);

//     fvPatchField<Type>::manipulateMatrix(eqn);
// }


template<class Type>
void immersedBoundaryFvPatchField<Type>::write(Ostream& os) const
{
    // to resolve the post-processing issues.  HJ, 1/Dec/2017
    fvPatchField<Type>::write(os);
    triValue_.writeEntry("triValue", os);
    triGrad_.writeEntry("triGradient", os);
    triValueFraction_.writeEntry("triValueFraction", os);
    os.writeKeyword("setDeadValue")
        << setDeadValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("deadValue")
        << deadValue_ << token::END_STATEMENT << nl;

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);

    // Write VTK on master only
    if (Pstream::master())
    {
        // Add parallel reduction of all faces and data to proc 0
        // and write the whola patch together
        
        // Write immersed boundary data as a vtk file
        autoPtr<surfaceWriter<Type> > writerPtr =
            surfaceWriter<Type>::New("vtk");

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
            surfaceWriterBase::FACE_DATA
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
