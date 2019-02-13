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

#include "immersedBoundaryFvcReconstruct.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "fv.H"
#include "surfaceInterpolate.H"
#include "fvcVolumeIntegrate.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcAverage.H"

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
reconstructIb
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

    // Note:
    // 1) Reconstruction is only available in cell centres: there is no need
    //    to invert the tensor on the boundary
    // 2) For boundaries, the only reconstructed data is the flux times
    //    normal.  Based on this guess, boundary conditions can adjust
    //    patch values
    // HJ, 12/Aug/2011

    GeometricField<GradType, fvPatchField, volMesh> fluxTimesNormal =
        surfaceSum((mesh.Sf()/mag(mesh.Sf()))*ssf);
      //  surfaceSum((mesh.Sf()/mesh.magSf())*ssf);

    // Note: hinv inverse must be used to stabilise the inverse on bad meshes
    // but it gives strange failures
    // HJ, 19/Aug/2015

    // Create the surface sum field
    GeometricField<SymmTensor<scalar>, fvPatchField, volMesh> surfaceSumSqrSf =
        surfaceSum(sqr(mesh.Sf())/mesh.magSf());
    Field<SymmTensor<scalar> > surfaceSumSqrSfIn = surfaceSumSqrSf.internalField();

    // Keep treck of the zero determinant cells
    dynamicLabelList openCells(20);

    // Check for zero determinant
    forAll(surfaceSumSqrSfIn, cellI)
    {
        if (mag(det(surfaceSumSqrSfIn[cellI])) < VSMALL)
        {
            // Set the field to a unit tensor
            surfaceSumSqrSfIn[cellI] = SymmTensor<scalar>::I;
            // Store the cell ID for later correction of the field
            openCells.append(cellI);
        }
    }

    reconField.internalField() =
    (
        inv
        (
           surfaceSumSqrSfIn
        )
      & fluxTimesNormal.internalField()
    );

    // Set the zero determinant cells reconstructed field to zero
    forAll(openCells, cellI)
    {
        reconField.internalField()[cellI] = GradType::zero;
    }

    reconField.boundaryField() = fluxTimesNormal.boundaryField();

    treconField().correctBoundaryConditions();

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
reconstructIb
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tvf
    (
        fvc::reconstructIb(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
