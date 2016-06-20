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

Description
    Specialisation of gaussDivScheme for vectors. Needed for implicit fvmDiv
    operator for block coupled systems.

\*---------------------------------------------------------------------------*/

#include "vectorGaussDivScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, scalar> > gaussDivScheme<vector>::fvmUDiv
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
    // Get weights
    surfaceScalarField weights = this->tinterpScheme_().weights(vf);

    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, scalar> > tbs
    (
        new BlockLduSystem<vector, scalar>(mesh)
    );
    BlockLduSystem<vector, scalar>& bs = tbs();
    scalarField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    // Internal field
    const vectorField& SfIn = mesh.Sf().internalField();
    const scalarField& wIn = weights.internalField();

    l = -wIn*SfIn;
    u = l + SfIn;
    bs.negSumDiag();

    // Boundary contributions
    forAll(vf.boundaryField(), patchI)
    {
        const fvPatchVectorField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& Sf = patch.Sf();
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
        const unallocLabelList& fc = patch.faceCells();

        const vectorField internalCoeffs(pf.valueInternalCoeffs(pw));

        // Diag contribution
        forAll(pf, faceI)
        {
            d[fc[faceI]] += cmptMultiply(internalCoeffs[faceI], Sf[faceI]);
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            const vectorField pcl = -pw*Sf;
            const vectorField pcu = pcl + Sf;

            // Coupling  contributions
            pcoupleLower -= pcl;
            pcoupleUpper -= pcu;
        }
        else
        {
            const vectorField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Boundary contribution
            forAll(pf, faceI)
            {
                source[fc[faceI]] -= boundaryCoeffs[faceI] & Sf[faceI];
            }
        }
    }

    // Interpolation schemes with corrections not accounted for

    return tbs;
}


template<>
tmp<BlockLduSystem<vector, scalar> > gaussDivScheme<vector>::fvmUDiv
(
    const surfaceScalarField& flux,
    const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
    // Get weights
    surfaceScalarField weights = this->tinterpScheme_().weights(vf);

    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, scalar> > tbs
    (
        new BlockLduSystem<vector, scalar>(mesh)
    );
    BlockLduSystem<vector, scalar>& bs = tbs();
    scalarField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    // Internal field
    const vectorField& SfIn = mesh.Sf().internalField();
    const scalarField& wIn = weights.internalField();
    const scalarField& fluxIn = flux.internalField();

    l = -wIn*fluxIn*SfIn;
    u = l + fluxIn*SfIn;
    bs.negSumDiag();

    // Boundary contributions
    forAll(vf.boundaryField(), patchI)
    {
        const fvPatchVectorField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& Sf = patch.Sf();
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
        const unallocLabelList& fc = patch.faceCells();

        const scalarField& pFlux = flux.boundaryField()[patchI];

        const vectorField internalCoeffs(pf.valueInternalCoeffs(pw));

        // Diag contribution
        forAll(pf, faceI)
        {
            d[fc[faceI]] += cmptMultiply
            (
                internalCoeffs[faceI],
                pFlux[faceI]*Sf[faceI]
            );
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            const vectorField pcl = -pw*pFlux*Sf;
            const vectorField pcu = pcl + pFlux*Sf;

            // Coupling  contributions
            pcoupleLower -= pcl;
            pcoupleUpper -= pcu;
        }
        else
        {
            const vectorField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Boundary contribution
            forAll(pf, faceI)
            {
                source[fc[faceI]] -=
                (
                    boundaryCoeffs[faceI] & (pFlux[faceI]*Sf[faceI])
                );
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
