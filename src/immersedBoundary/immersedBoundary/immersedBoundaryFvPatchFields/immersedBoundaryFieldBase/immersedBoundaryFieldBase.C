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

#include "immersedBoundaryFieldBase.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedBoundaryFieldBase<Type>::immersedBoundaryFieldBase
(
    const fvPatch& p,
    const bool setDeadValue,
    const Type deadValue
)
:
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    setDeadValue_(setDeadValue),
    deadValue_(deadValue)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedBoundaryFieldBase<Type>::setDeadValues
(
    Field<Type>& psiI
) const
{
    // Fix the value in dead cells
    if (setDeadValue_)
    {
        const labelList& dc = ibPatch_.ibPolyPatch().deadCells();

        forAll (dc, dcI)
        {
            psiI[dc[dcI]] = deadValue_;
        }
    }
}


template<class Type>
void Foam::immersedBoundaryFieldBase<Type>::setDeadValues
(
    fvMatrix<Type>& matrix
) const
{
    // Fix the value in dead cells
    if (setDeadValue_)
    {
        const labelList& dc = ibPatch_.ibPolyPatch().deadCells();

        // Boost the diagonal of dead cells by the volume ratio
        // Volume ratio is set to SMALL; revert for diagonal
        // This should also guarantee strong diagonal dominance.
        // HJ, 19/Jun/2019

        scalarField& diag = matrix.diag();

        forAll (dc, dcI)
        {
            diag[dc[dcI]] *= GREAT;
        }

        // Set values
        matrix.setValues
        (
            dc,
            Field<Type>(dc.size(), deadValue_)
        );
    }
}


template<class Type>
void Foam::immersedBoundaryFieldBase<Type>::writeDeadData(Ostream& os) const
{
    os.writeKeyword("setDeadValue")
        << setDeadValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("deadValue")
        << deadValue_ << token::END_STATEMENT << nl;
}


template<class Type>
void Foam::immersedBoundaryFieldBase<Type>::writeField
(
    const fvPatchField<Type>& f
) const
{
    // Write immersed boundary data as a vtk file
    autoPtr<surfaceWriter> writerPtr = surfaceWriter::New("vtk");

    // Get the intersected patch
    const standAlonePatch& ts = ibPatch_.ibPolyPatch().ibPatch();

    if (Pstream::parRun())
    {
        // Gather points, fields and faces

        // Points
        List<pointField> procPoints(Pstream::nProcs());
        procPoints[Pstream::myProcNo()] = ts.points();
        Pstream::gatherList(procPoints);

        // Fields
        List<Field<Type> > procFields(Pstream::nProcs());
        procFields[Pstream::myProcNo()] = f;
        Pstream::gatherList(procFields);

        // Faces
        List<faceList> procFaces(Pstream::nProcs());
        procFaces[Pstream::myProcNo()] = ts;
        Pstream::gatherList(procFaces);

        if (Pstream::master())
        {
            // Assemble unique lists to correspond to a single surface

            // Count points and faces; currently unmerged
            label nAllPoints = 0;
            label nAllFaces = 0;

            forAll (procPoints, procI)
            {
                nAllPoints += procPoints[procI].size();
                nAllFaces += procFaces[procI].size();
            }

            pointField allPoints(nAllPoints);
            faceList allFaces(nAllFaces);
            Field<Type> completeField(nAllFaces);

            // Reset counters
            nAllPoints = 0;
            nAllFaces = 0;

            forAll (procPoints, procI)
            {
                const pointField& curPoints = procPoints[procI];
                labelList renumberPoints(curPoints.size());

                forAll (curPoints, cpI)
                {
                    allPoints[nAllPoints] = curPoints[cpI];
                    renumberPoints[cpI] = nAllPoints;
                    nAllPoints++;
                }

                const faceList& curFaces = procFaces[procI];
                const Field<Type>& curField = procFields[procI];

                forAll (curFaces, cfI)
                {
                    // Point labels in faces need to be renumbered with respect
                    // to the size of the size of the previous processore patch
                    const face& oldFace = curFaces[cfI];

                    // Make a copy of face to renumber
                    face renumberedFace(oldFace.size());

                    forAll (oldFace, fpI)
                    {
                        renumberedFace[fpI] = renumberPoints[oldFace[fpI]];
                    }

                    allFaces[nAllFaces] = renumberedFace;
                    completeField[nAllFaces] = curField[cfI];
                    nAllFaces++;
                }
            }

            if (nAllPoints != allPoints.size() || nAllFaces != allFaces.size())
            {
                FatalErrorInFunction
                    << "Problem with merge of immersed boundary patch data"
                    << abort(FatalError);
            }

            writerPtr->write
            (
                f.dimensionedInternalField().path(),
                ibPatch_.name(),
                allPoints,
                allFaces,
                f.dimensionedInternalField().name(),
                completeField,
                false, // FACE_DATA
                false  // verbose
            );
        }
    }
    else
    {
        writerPtr->write
        (
            f.dimensionedInternalField().path(),
            ibPatch_.name(),
            ts.points(),
            ts,
            f.dimensionedInternalField().name(),
            f,
            false, // FACE_DATA
            false  // verbose
        );
    }
}


// ************************************************************************* //
