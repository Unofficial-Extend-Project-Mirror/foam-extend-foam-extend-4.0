/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "gaussDivScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
>
gaussDivScheme<Type>::fvcDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
    > tDiv
    (
        fvc::surfaceIntegrate
        (
            this->mesh_.Sf() & this->tinterpScheme_().interpolate(vf)
        )
    );

    tDiv().rename("div(" + vf.name() + ')');

    return tDiv;
}


template<class Type>
tmp<blockVectorMatrix> gaussDivScheme<Type>::fvmDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(vf);
    const scalarField& wIn = tweights().internalField();

    const fvMesh& mesh = vf.mesh();

    tmp<blockVectorMatrix> tbm
    (
        new blockVectorMatrix
        (
           mesh
        )
    );
    blockVectorMatrix& bm = tbm();

    // Grab ldu parts of block matrix as linear always
    typename CoeffField<vector>::linearTypeField& u = bm.upper().asLinear();
    typename CoeffField<vector>::linearTypeField& l = bm.lower().asLinear();

    const vectorField& SfIn = mesh.Sf().internalField();

    l = -wIn*SfIn;
    u = l + SfIn;
    bm.negSumDiag();

    // Interpolation schemes with corrections not accounted for

    return tbm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
