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

#include "explicitAdjConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type> >
explicitAdjConvectionScheme<Type>::fvmAdjDiv
(
    const volVectorField& Up,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*Up.dimensions()*vf.dimensions()/dimLength
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    // Matrix consistency
    fvm.diag() = 0;

    fvm += this->fvcAdjDiv(Up, vf);

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
explicitAdjConvectionScheme<Type>::fvcAdjDiv
(
    const volVectorField& Up,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Return Up & grad(vf).T()
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAdjConvection
    (
        fvc::grad(vf) & Up
    );

    tAdjConvection().rename
    (
        "adjConvection(" + Up.name() + ',' + vf.name() + ')'
    );

    return tAdjConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
