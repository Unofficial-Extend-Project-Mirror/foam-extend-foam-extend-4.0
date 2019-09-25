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

#include "fvcReconstruct.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
reconstruct
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate(" + ssf.name() + ')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            ssf.dimensions()/dimArea,
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& reconField =
        treconField();

    // Notes regarding boundary:
    // 1. Reconstruction is only available in cell centres: there is no need
    //    to invert the tensor on the boundary
    // 2. For boundaries, the only reconstructed data is the flux times
    //    normal. Based on this guess, boundary conditions can adjust
    //    patch values. HJ, 12/Aug/2011

    // Notes regarding procedure:
    // 1. hinv inverse must be used to stabilise the inverse on bad meshes
    //    but it gives strange failures because it unnecessarily avoids
    //    performing ordinary inverse for meshes with reasonably sized
    //    determinant (e.g. if SfSf/magSf is small). HJ, 19/Aug/2015
    // 2. hinv has been stabilised now. HJ, 22/Mar/2019
    // 3. But we still need to make sure that the determinant is not extremely
    //    small, which may happen for extremely small meshes. We avoid this
    //    issue by dividing the reconstruction equation with magSf^2 (instead of
    //    magSf), which basically makes the dyadic tensor that we need to invert
    //    dimensionless. VV, 13/Jun/2019

    // Here's a short derivation in a Latex--like notation, where:
    // - Sf is the surface area vector
    // - uf is the face velocity factor (or field to be reconstructed)
    // - Ff is the face flux
    // - magSf is the magnitude of the surface area vector
    // - uP is the velocity field in the cell centre
    // - G is the surface area dyadic tensor
    // - Fn is the vector representing surface sum of directional fluxes
    // - \dprod is a dot product

    // 1. Sf \dprod uf = Ff
    //    Multiply Eq (1) with Sf/magSf^2
    // 2. \frac{Sf Sf}{magSf^2} \dprod uf = \frac{Sf Ff}{magSf^2}
    //    Sum Eq (2) over all the faces
    // 3. \sum_f(\frac{Sf Sf}{magSf^2} \dprod uf) = \sum_f(\frac{Sf Ff}{magSf^2})
    //    Assume first order extrapolation of uf, e.g. uP = uf
    // 4. \sum_f(\frac{Sf Sf}{magSf^2}) \dprod uP) = \sum_f(\frac{Sf Ff}{magSf^2})
    //    Use shorthand notation
    // 5. G \dprod uP = Fn
    // 6. uP = G^-1 \dprod Fn

    // Fn -> fluxTimesNormal
    // G  -> G

    // Calculate sum of the directional fluxes
    const surfaceScalarField magSfSqr =
        sqr(mesh.magSf() + dimensionedScalar("vsmall", dimArea, 1e-100));

    const GeometricField<GradType, fvPatchField, volMesh> fluxTimesNormal =
        surfaceSum((mesh.Sf()/magSfSqr)*ssf);

    // Calculate the G tensor
    const volSymmTensorField G = surfaceSum(sqr(mesh.Sf())/magSfSqr);

    // Finally calculate the reconstructed field using hinv for stabilisation on
    // really bad fvMesh bits (uses ordinary inverse most of the time, see
    // tensor.C)
    reconField.internalField() =
        hinv(G.internalField()) & fluxTimesNormal.internalField();

    // Boundary value update
    reconField.boundaryField() = fluxTimesNormal.boundaryField();
    reconField.correctBoundaryConditions();

    return treconField;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
reconstruct
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tvf
    (
        fvc::reconstruct(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
