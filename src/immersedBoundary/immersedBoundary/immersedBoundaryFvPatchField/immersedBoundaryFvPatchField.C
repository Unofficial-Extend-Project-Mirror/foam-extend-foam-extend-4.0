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
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void immersedBoundaryFvPatchField<Type>::updateIbValues() const
{
    // Interpolate the values form tri surface using nearest triangle
    ibValue_ = Field<Type>(refValue, ibPatch_.nearestTri());
    ibGrad_ = Field<Type>(refGrad, ibPatch_.nearestTri());
    ibValueFraction_ = scalarField(refValueFraction, ibPatch_.nearestTri());
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::imposeDeadCondition()
{
    const labelList& dc = ibPatch_.deadCells();

    // Get non-const access to internal field
    Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

    forAll (dc, dcI)
    {
        psiI[dc[dcI]] = deadCellValue_;
    }
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::correctDiag
(
    fvMatrix<Type>& eqn
) const
{
    scalarField& Diag = eqn.diag();

    const labelList& dce = ibPatch_.deadCellsExt();

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    if (dce.size() < Diag.size())
    {
        liveDiag = gSumMag(Diag)/(Diag.size() - dce.size());

        // Correct for sign
        liveDiag *= sign(gMax(Diag));
    }

    forAll (dce, cellI)
    {
        if (mag(Diag[dce[cellI]]) < SMALL)
        {
            Diag[dce[cellI]] = liveDiag;
        }
    }
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::correctOffDiag
(
    fvMatrix<Type>& eqn
) const
{
    // Calculate gradient contribution
    const labelList& ibFaces = ibPatch_.ibFaces();
    const labelList& ibFaceCells = ibPatch_.ibFaceCells();

    const scalarField& ibGamma = ibPatch_.gamma().internalField();

    const unallocLabelList& own = mesh_.owner();
    const unallocLabelList& nei = mesh_.neighbour();

    // Get delta coefficients
    const surfaceScalarField& dc = mesh_.deltaCoeffs();
    const scalarField& dcI = dc.internalField();

    if (eqn.symmetric())
    {
        scalarField& diag = eqn.diag();
        scalarField& upper = eqn.upper();
        Field<Type>& source = eqn.source();

//         Info<< "Symmetric correctOffDiag for field "
//             << this->dimensionedInternalField().name() << endl;

        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {
                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];

                    source[own[curFace]] +=
                        upper[curFace]*ibGrad_[ibFaceCells[faceI]]
                        /dcI[curFace];
                }
                else
                {
                    diag[nei[curFace]] += upper[curFace];

                    source[nei[curFace]] -=
                        upper[curFace]*ibGrad_[ibFaceCells[faceI]]
                        /dcI[curFace];
                }

                upper[curFace] = 0;
            }
            else
            {
                // else MPH
                label patchi = mesh_.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh_.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    // Check if the live cell is on local or neighbour side
                    // HJ, 7/Dec/2012
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        {
                            source[ibFaceCells[faceI]] +=
                                ibGrad_[ibFaceCells[faceI]]
                                /dc.boundaryField()[patchi][patchFacei];
                        }
                    }
                }
            }
        }
    }
    else if (eqn.asymmetric())
    {
        scalarField& diag = eqn.diag();
        scalarField& upper = eqn.upper();
        scalarField& lower = eqn.lower();
        Field<Type>& source = eqn.source();

//         Info<< "Asymmetric correctOffDiag for field "
//             << this->dimensionedInternalField().name() << endl;

        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {
                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];

                    source[own[curFace]] +=
                        upper[curFace]*ibGrad_[ibFaceCells[faceI]]/dcI[faceI];
                }
                else
                {
                    diag[nei[curFace]] += lower[curFace];

                    source[nei[curFace]] -=
                        lower[curFace]*ibGrad_[ibFaceCells[faceI]]/dcI[faceI];
                }

                upper[curFace] = 0;
                lower[curFace] = 0;
            }
            else
            {
                // else MPH
                label patchi = mesh_.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh_.boundaryMesh()[patchi].whichFace(curFace);

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    // Check if the live cell is on local or neighbour side
                    // HJ, 7/Dec/2012
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        {
                            source[ibFaceCells[faceI]] +=
                                ibGrad_[ibFaceCells[faceI]]
                                /dc.boundaryField()[patchi][patchFacei];
                        }
                    }
                }
            }
        }
    }

    // Note: potentially deal with face flux correction ptr.
    // HJ, 16/Apr/2012
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool immersedBoundaryFvPatchField<Type>::motionUpdateRequired() const
{
    if
    (
        ibPatch_.movingIb()
     || ibPatch_.boundaryMesh().mesh().moving()
    )
    {
        if
        (
            ibUpdateTimeIndex_
         != ibPatch_.boundaryMesh().mesh().time().timeIndex()
        )
        {
            // Mesh is moving and current time has not been updated
            return true;
        }
    }

    return false;
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::motionUpdate() const
{
    // Motion update, clear data related to immersed boundary points
    ibValue_.clear();
    ibGrad_.clear();

    // Record motion update time
    ibUpdateTimeIndex_ = ibPatch_.boundaryMesh().mesh().time().timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    refGrad_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    refValueFraction_(false),
    setDeadCellValue_(false),
    deadCellValue_(pTraits<Type>::zero),
    ibValue_(),
    ibGrad_()
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_("refValue", dict, ibPatch_.ibMesh().size()),
    refGrad_("refGradient", dict, ibPatch_.ibMesh().size()),
    refValueFraction_(dict.lookup("valueFraction")),
    setDeadCellValue_(dict.lookup("setDeadCellValue")),
    deadCellValue_(pTraits<Type>(dict.lookup("deadCellValue"))),
    ibValue_(),
    ibGrad_()
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
    fvPatchField<Type>(ptf, p, iF, mapper),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    refValueFraction_(ptf.refValueFraction()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibValue_(),
    ibGrad_()
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
}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    ibPatch_(ptf.ibPatch()),
    mesh_(ptf.patch().boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    refValueFraction_(ptf.refValueFraction()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibValue_(),
    ibGrad_()
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    ibPatch_(ptf.ibPatch()),
    mesh_(ptf.patch().boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    refValueFraction_(ptf.refValueFraction()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibValue_(),
    ibGrad_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Field<Type>& immersedBoundaryFvPatchField<Type>::ibValue() const
{
    // Note: on a moving mesh, the intersection has changed and
    // ibValue and ibGrad fields should be cleared and recalculated.
    // HJ, 17/Oct/2012
    if (this->motionUpdateRequired())
    {
        motionUpdate();
    }

    if (ibValue_.empty())
    {
        this->updateIbValues();
    }

    return ibValue_;
}


template<class Type>
const Field<Type>& immersedBoundaryFvPatchField<Type>::ibGrad() const
{
    // Note: on a moving mesh, the intersection has changed and
    // ibValue and ibGrad fields should be cleared and recalculated.
    // HJ, 17/Oct/2012
    if (this->motionUpdateRequired())
    {
        motionUpdate();
    }

    if (ibGrad_.empty())
    {
        this->updateIbValues();
    }

    return ibGrad_;
}


template<class Type>
tmp<Field<Type> > immersedBoundaryFvPatchField<Type>::ibCellValue() const
{
    // Collect IB cell values
    tmp<Field<Type> > tibcv
    (
        new Field<Type>
        (
            this->internalField(),
            ibPatch_.ibCells()
        )
    );

    return tibcv;
}


template<class Type>
tmp<Field<Type> >
immersedBoundaryFvPatchField<Type>::ibSamplingPointValue() const
{
    return ibPatch_.toSamplingPoints(this->internalField());
}


template<class Type>
tmp<Field<Type> > immersedBoundaryFvPatchField<Type>::triValue() const
{
    return ibPatch_.toTriFaces(this->ibValue());
}


template<class Type>
tmp<Field<Type> > immersedBoundaryFvPatchField<Type>::triGrad() const
{
    return ibPatch_.toTriFaces(this->ibGrad());
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::updateCoeffs()
{
    if (this->fixesValue())
    {
        this->imposeDirichletCondition();
    }
    else
    {
        this->imposeNeumannCondition();
    }

    // Fix the value in dead cells
    if (setDeadCellValue_)
    {
        this->imposeDeadCondition();
    }

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes
)
{
    this->updateCoeffs();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    // Note
    // Since the boundary condition is performed by data fitting with the
    // internal field, fitting must be performed both on updateCoeffs
    // and on evaluate (internal field has changed in the meantime).
    // Bug fix, Zeljko Tukovic, 21/Jun/2012
//     this->updateCoeffs();

    fvPatchField<Type>::evaluate();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& eqn
)
{
    this->initEvaluate();

    // Build matrix diagonal for cells where it is missing
    this->correctDiag(eqn);

    // For Neumann boundary condition, manipulate matrix off-diagonal
    if (!this->fixesValue())
    {
        this->correctOffDiag(eqn);
    }

    // Set values in IB cells
    Field<Type> polyPsi(eqn.psi(), ibPatch_.ibCells());
    eqn.setValues(ibPatch_.ibCells(), polyPsi);

    // Correct equation for dead cells
    Field<Type> deadCellsPsi
    (
        ibPatch_.deadCells().size(),
        deadCellValue_
    );
    eqn.setValues(ibPatch_.deadCells(), deadCellsPsi);

    fvPatchField<Type>::manipulateMatrix(eqn);
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    refValueFraction_.writeEntry("valueFraction", os);
    os.writeKeyword("setDeadCellValue")
        << setDeadCellValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("deadCellValue")
        << deadCellValue_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);

    // Write VTK on master only
    // HJ. fix here!!!

    // if (Pstream::master())
    // {
    //     // Write immersed boundary data as a vtk file
    //     autoPtr<surfaceWriter<Type> > writerPtr =
    //         surfaceWriter<Type>::New("vtk");

    //     const triSurface& ts = ibPatch_.ibMesh();

    //     // Make a face list for writing
    //     faceList f(ts.size());
    //     forAll (ts, faceI)
    //     {
    //         f[faceI] = ts[faceI].triFaceFace();
    //     }

    //     writerPtr->write
    //     (
    //         this->dimensionedInternalField().path(),
    //         ibPatch_.name(),
    //         ts.points(),
    //         f,
    //         this->dimensionedInternalField().name(),
    //         triData,
    //         surfaceWriterBase::FACE_DATA
    //     );
    // }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
