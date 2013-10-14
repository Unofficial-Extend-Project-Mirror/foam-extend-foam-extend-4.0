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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/


#include "TractionPointPatchVectorField.H"
#include "constraints.H"
#include "tetFemMatrix.H"
#include "primitivePatch.H"
#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
void TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>
::checkFieldSize() const
{
    if
    (
        this->size() != this->patch().size()
     || traction_.size() != this->patch().nFaces()
     || pressure_.size() != this->patch().nFaces()
    )
    {
        FatalErrorIn
        (
            "void TractionPointPatchVectorField::checkField() const"
        )   << "field does not correspond to patch. " << endl
            << "Field size: " << this->size()
            << " traction size: " << traction_.size()
            << " pressure size: " << pressure_.size()
            << " patch size: " << this->patch().size()
            << " patch number of faces: " << this->patch().nFaces()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const DimensionedField<vector, tetPointMesh>& iF
)
:
    PatchField<vector>(p, iF),
    traction_(this->patch().nFaces()),
    pressure_(this->patch().nFaces())
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const DimensionedField<vector, tetPointMesh>& iF,
    const vectorField& v,
    const scalarField& vf
)
:
    PatchField<vector>(p, iF),
    traction_(v),
    pressure_(vf)
{
    checkFieldSize();
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
TractionPointPatchVectorField
(
    const PointPatch& p,
    const DimensionedField<vector, tetPointMesh>& iF,
    const dictionary& dict
)
:
    PatchField<vector>(p, iF),
    traction_("traction", dict, this->patch().nFaces()),
    pressure_("pressure", dict, this->patch().nFaces())
{
    this->updateBoundaryField();
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
TractionPointPatchVectorField
(
    const TractionPointPatchVectorField
        <PatchField, Mesh, PointPatch, MatrixType>& ptf,
    const PointPatch& p,
    const DimensionedField<vector, tetPointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    PatchField<vector>(p, iF),
    traction_(ptf.traction_, mapper),
    pressure_(ptf.pressure_, mapper)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
TractionPointPatchVectorField
(
    const TractionPointPatchVectorField
    <PatchField, Mesh, PointPatch, MatrixType>& ptf,
    const DimensionedField<vector, tetPointMesh>& iF
)
:
    PatchField<vector>(ptf, iF),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Grab the values using rmap
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
void
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::rmap
(
    const PointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, vector>& ptf,
    const labelList& addr
)
{
    const TractionPointPatchVectorField& mptf =
        refCast<const TractionPointPatchVectorField>(ptf);

    traction_.rmap(mptf.traction_, addr);
    pressure_.rmap(mptf.pressure_, addr);
}


// Set boundary condition to matrix
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
void TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
addBoundarySourceDiag
(
    MatrixType<vector>& matrix
) const
{
    const vectorField& points = this->patch().localPoints();

    const labelList& meshPoints = this->patch().meshPoints();
    const vectorField& pointNormals = this->patch().pointNormals();
    const label nFaces = this->patch().nFaces();

    vectorField& source = matrix.source();

    for (label faceI = 0; faceI < nFaces; faceI++)
    {
        List<triFace> faceTriangles =
            this->patch().faceTriangles(faceI);

        forAll (faceTriangles, triI)
        {
            const triFace& tri = faceTriangles[triI];

            source[meshPoints[tri[0]]] +=
                tri.mag(points)*
                (
                    traction_[faceI]/6.0
                  - pressure_[faceI]*pointNormals[tri[0]]/6.0
                  + traction_[faceI]/12.0
                  - pressure_[faceI]*pointNormals[tri[1]]/12.0
                  + traction_[faceI]/12.0
                  - pressure_[faceI]*pointNormals[tri[2]]/12.0
                );

            source[meshPoints[tri[1]]] +=
                tri.mag(points)*
                (
                    traction_[faceI]/6.0
                  - pressure_[faceI]*pointNormals[tri[1]]/6.0
                  + traction_[faceI]/12.0
                  - pressure_[faceI]*pointNormals[tri[2]]/12.0
                  + traction_[faceI]/12.0
                  - pressure_[faceI]*pointNormals[tri[0]]/12.0
                );

            source[meshPoints[tri[2]]] +=
                tri.mag(points)*
                (
                    traction_[faceI]/6.0
                  - pressure_[faceI]*pointNormals[tri[2]]/6.0
                  + traction_[faceI]/12.0
                  - pressure_[faceI]*pointNormals[tri[0]]/12.0
                  + traction_[faceI]/12.0
                  - pressure_[faceI]*pointNormals[tri[1]]/12.0
                );
        }
    }
}


// Write
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType
>
void
TractionPointPatchVectorField<PatchField, Mesh, PointPatch, MatrixType>::
write(Ostream& os) const
{
    PatchField<vector>::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
