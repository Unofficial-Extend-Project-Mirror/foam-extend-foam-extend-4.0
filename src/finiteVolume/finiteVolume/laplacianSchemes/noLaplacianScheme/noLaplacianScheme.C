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
    
\*---------------------------------------------------------------------------*/

#include "noLaplacianScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
noLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const surfaceScalarField& deltaCoeffs = this->mesh().deltaCoeffs();
    const surfaceScalarField& magSf = this->mesh().magSf();

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "laplacian("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            vf.mesh(),
            dimensioned<Type>
            (
                "0",
                deltaCoeffs.dimensions()*magSf.dimensions()*vf.dimensions(),
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
noLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const surfaceScalarField& deltaCoeffs = this->mesh().deltaCoeffs();
    const surfaceScalarField& magSf = this->mesh().magSf();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gamma.dimensions()*
            magSf.dimensions()*vf.dimensions()
        )
    );

    return tfvm;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
noLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const surfaceScalarField& deltaCoeffs = this->mesh().deltaCoeffs();
    const surfaceScalarField& magSf = this->mesh().magSf();

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "laplacian("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            vf.mesh(),
            dimensioned<Type>
            (
                "0",
                deltaCoeffs.dimensions()*gamma.dimensions()*
                magSf.dimensions()*vf.dimensions(),
                pTraits<Type>::zero
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
