/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "fvmDiv.H"
#include "fvMesh.H"
#include "convectionScheme.H"
#include "divScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type> >
div
(
    const surfaceScalarField& flux,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().schemesDict().divScheme(name)
    )().fvmDiv(flux, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
div
(
    const tmp<surfaceScalarField>& tflux,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type> > Div(fvm::div(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<fvMatrix<Type> >
div
(
    const surfaceScalarField& flux,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::div(flux, vf, "div("+flux.name()+','+vf.name()+')');
}


template<class Type>
tmp<fvMatrix<Type> >
div
(
    const tmp<surfaceScalarField>& tflux,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > Div(fvm::div(tflux(), vf));
    tflux.clear();
    return Div;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> UDiv
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().divScheme(name)
    )().fvmUDiv(vf);
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> UDiv
(
    const surfaceScalarField& flux,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().divScheme(name)
    )().fvmUDiv(flux, vf);
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> UDiv
(
    const tmp<surfaceScalarField>& tflux,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    >
    Div(fvm::UDiv(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> UDiv
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::UDiv
    (
        vf,
        "div(" + vf.name() + ')'
    );
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> UDiv
(
    const surfaceScalarField& flux,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::UDiv
    (
        flux,
        vf,
        "div(" + vf.name() + ')'
    );
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> UDiv
(
    const tmp<surfaceScalarField>& tflux,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp
    <
        BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
    >
    Div(fvm::UDiv(tflux(), vf));
    tflux.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
