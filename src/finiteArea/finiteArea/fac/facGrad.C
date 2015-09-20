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

#include "facGrad.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "facEdgeIntegrate.H"
#include "faMesh.H"
#include "faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
grad
(
    const GeometricField<Type, faePatchField, edgeMesh>& ssf
)
{
    return fac::edgeIntegrate(ssf.mesh().Sf() * ssf);
}

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, faPatchField, areaMesh
    >
>
grad
(
    const tmp<GeometricField<Type, faePatchField, edgeMesh> >& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, faPatchField, areaMesh> > Grad
    (
        fac::grad(tssf())
    );
    tssf.clear();
    return Grad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, faPatchField, areaMesh
    >
>
grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::gradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().gradScheme(name)
    )().grad(vf);
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, faPatchField, areaMesh
    >
>
grad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type, faPatchField, areaMesh
        >
    > tGrad
    (
        fac::grad(tvf(), name)
    );
    tvf.clear();
    return tGrad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, faPatchField, areaMesh
    >
>
grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::grad(vf, "grad(" + vf.name() + ')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, faPatchField, areaMesh
    >
>
grad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, faPatchField, areaMesh> > Grad
    (
        fac::grad(tvf())
    );
    tvf.clear();
    return Grad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
