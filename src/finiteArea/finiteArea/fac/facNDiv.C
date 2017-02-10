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

#include "facNDiv.H"
#include "faMesh.H"
#include "facEdgeIntegrate.H"
#include "faDivScheme.H"
#include "faConvectionScheme.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const GeometricField<Type, faePatchField, edgeMesh>& ssf
)
{
    const areaVectorField& n = ssf.mesh().faceAreaNormals();

    tmp<GeometricField<Type, faPatchField, areaMesh> > v =
        fac::edgeIntegrate(ssf);

    v.internalField() = n*(n & v.internalField());
    v.correctBoundaryConditions();

    return v;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const tmp<GeometricField<Type, faePatchField, edgeMesh> >& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div(fac::ndiv(tssf()));
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
ndiv
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    const areaVectorField& n = vf.mesh().faceAreaNormals();

    tmp<GeometricField<Type, faPatchField, areaMesh> > tDiv
    (
        fa::divScheme<Type>::New
            (
                vf.mesh(), vf.mesh().schemesDict().divScheme(name)
            )().facDiv(vf)
    );
    GeometricField<Type, faPatchField, areaMesh>& Div = tDiv();

    Div.internalField() = n*(n & Div.internalField());
    Div.correctBoundaryConditions();

    return tDiv;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
ndiv
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh> > Div
    (
        fac::ndiv(tvvf(), name)
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
ndiv
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::ndiv(vf, "div("+vf.name()+')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
ndiv
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh> > Div
    (
        fac::ndiv(tvvf())
    );

    tvvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    const areaVectorField& n = vf.mesh().faceAreaNormals();

    tmp<GeometricField<Type, faPatchField, areaMesh> > tDiv
    (
        fa::convectionScheme<Type>::New
        (
            vf.mesh(),
            flux,
            vf.mesh().schemesDict().divScheme(name)
        )().facDiv(flux, vf)
    );

    GeometricField<Type, faPatchField, areaMesh>& Div = tDiv();

    Div.internalField() = n*(n &Div.internalField());
    Div.correctBoundaryConditions();

    return tDiv;

}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::ndiv(tflux(), vf, name)
    );
    tflux.clear();

    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::ndiv(flux, tvf(), name)
    );
    tvf.clear();

    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::ndiv(tflux(), tvf(), name)
    );
    tflux.clear();
    tvf.clear();

    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::ndiv
    (
        flux, vf, "div("+flux.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::ndiv(tflux(), vf)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::ndiv(flux, tvf())
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ndiv
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Div
    (
        fac::ndiv(tflux(), tvf())
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
