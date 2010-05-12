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

Description
    

\*---------------------------------------------------------------------------*/

#include "facDiv.H"
#include "faMesh.H"
#include "facEdgeIntegrate.H"
#include "faDivScheme.H"
#include "faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const GeometricField<Type, faePatchField, edgeMesh>& ssf
)
{
    return fac::edgeIntegrate(ssf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const tmp<GeometricField<Type, faePatchField, edgeMesh> >& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div(fac::div(tssf()));
    tssf.clear();
    return Div;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::divScheme<Type>::New
    (
        vf.mesh(), vf.mesh().divScheme(name)
    )().facDiv(vf);
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh> > Div
    (
        fac::div(tvvf(), name)
    );
    tvvf.clear();
    return Div;
}

template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh> > Div(fac::div(tvvf()));
    tvvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().divScheme(name)
    )().facDiv(flux, vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::div(tflux(), vf, name)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::div(flux, tvf(), name)
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::div(tflux(), tvf(), name)
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::div
    (
        flux, vf, "div("+flux.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::div(tflux(), vf)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::div(flux, tvf())
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
div
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::div(tflux(), tvf())
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
