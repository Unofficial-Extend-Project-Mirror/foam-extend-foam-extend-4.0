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

Description
    Specialisation of gaussGrad for scalars. Needed for implicit fvmGrad
    operator for block coupled systems.

\*---------------------------------------------------------------------------*/

#include "scalarGaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, vector> > gaussGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(vf);
    const scalarField& wIn = tweights().internalField();

    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(mesh)
    );
    BlockLduSystem<vector, vector>& bs = tbs();
    vectorField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    typename CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    typename CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    typename CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    const vectorField& SfIn = mesh.Sf().internalField();

    l = -wIn*SfIn;
    u = l + SfIn;
    bs.negSumDiag();

    // Boundary contributions
    forAll(vf.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& Sf = patch.Sf();
        const fvsPatchScalarField& pw = tweights().boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));

        // Diag contribution
        forAll(pf, faceI)
        {
            d[fc[faceI]] += internalCoeffs[faceI]*Sf[faceI];
        }

        if (patch.coupled())
        {
            typename CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            typename CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            const vectorField pcl = -pw*Sf;
            const vectorField pcu = pcl + Sf;

            // Coupling  contributions
            pcoupleLower -= pcl;
            pcoupleUpper -= pcu;
        }
        else
        {
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Boundary contribution
            forAll(pf, faceI)
            {
                source[fc[faceI]] -= boundaryCoeffs[faceI]*Sf[faceI];
            }
        }
    }

    // Interpolation schemes with corrections not accounted for

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
