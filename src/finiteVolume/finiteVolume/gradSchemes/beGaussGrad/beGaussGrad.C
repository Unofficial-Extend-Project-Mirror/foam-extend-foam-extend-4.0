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

#include "beGaussGrad.H"
#include "zeroGradientFvPatchField.H"

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
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
beGaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tssf =
        tinterpScheme_().interpolate(vsf);

    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf();

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    // Calculate the gradient on the internal field
    forAll(owner, facei)
    {
        GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

    // Assemble boundary contribution
    // If a patch is coupled or fixes value, regular gradient calculation
    // is used.  For all other patch types, one-sided implicit extrapolation
    // is applied instead.  HJ, 18/Sep/2006 Experimental!

    tensorField bct(mesh.V().field()*tensor(sphericalTensor::I));
    const vectorField& Cin = mesh.C();
    const Field<Type>& vsfIn = vsf.internalField();

    forAll (mesh.boundary(), patchi)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const vectorField& pfC = mesh.C().boundaryField()[patchi];

        if
        (
            vsf.boundaryField()[patchi].coupled()
         || vsf.boundaryField()[patchi].fixesValue()
        )
        {
            // Regular calculation, using patch value
            const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

            forAll(mesh.boundary()[patchi], facei)
            {
                igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
            }
        }
        else
        {
            // One-sided implicit extrapolation, using cell value
            forAll(mesh.boundary()[patchi], facei)
            {
                igGrad[pFaceCells[facei]] +=
                    pSf[facei]*vsfIn[pFaceCells[facei]];

                // Accumulate face dot products for inversion
                bct[pFaceCells[facei]] -=
                    pSf[facei]*(pfC[facei] - Cin[pFaceCells[facei]]);
            }
        }
    }

    // Calculate corrected gradient and evaluate boundary values
    // Note: boundary correction for forced snGrad is removed
    // HJ, 19/Sep/2006
    igGrad = (igGrad & inv(bct));
    gGrad.correctBoundaryConditions();

    return tgGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
