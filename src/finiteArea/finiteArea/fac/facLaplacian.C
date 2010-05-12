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
    Namespace of functions to calculate explicit derivatives.
    Time derivatives are calculated using Euler-implicit, backward differencing
    or Crank-Nicholson. Spatial derivatives are calculated using Gauss' Theorem.

\*---------------------------------------------------------------------------*/

#include "facLaplacian.H"
#include "faMesh.H"
#include "faLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().facLaplacian(vf);
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian(vf, "laplacian(" + vf.name() + ')');
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const dimensionedScalar& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return gamma*fac::laplacian(vf, name);
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const dimensionedScalar& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const dimensionedScalar& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return gamma*fac::laplacian
    (
        vf, "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const dimensionedScalar& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const areaScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().facLaplacian(gamma, vf);
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<areaScalarField>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const areaScalarField& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<areaScalarField>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const areaScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<areaScalarField>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian
    (
        tgamma,
        vf,
        "laplacian(" + tgamma().name() + ',' + vf.name() + ')'
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const areaScalarField& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    return fac::laplacian
    (
        gamma,
        tvf,
        "laplacian(" + gamma.name() + ',' + tvf().name() + ')'
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<areaScalarField>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    return fac::laplacian
    (
        tgamma,
        tvf,
        "laplacian(" + tgamma().name() + ',' + tvf().name() + ')'
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const edgeScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().facLaplacian(gamma, vf);
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<edgeScalarField>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const edgeScalarField& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
(
    const tmp<edgeScalarField>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const edgeScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<edgeScalarField>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const edgeScalarField& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
(
    const tmp<edgeScalarField>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


/*
template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const areaTensorField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = tvf().mesh();
    return fac::laplacian
    (
        (
            mesh.Sf() & fac::interpolate(tgamma) & mesh.Sf()
        )/sqr(mesh.magSf()),
        tvf
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const areaTensorField& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<areaTensorField>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<areaTensorField>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
(
    const edgeTensorField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = tvf().mesh();

    return fac::laplacian
    (
        (
            mesh.Sf() &
            gamma
          & mesh.Sf()
        )/sqr(mesh.magSf()),
        vf
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const edgeTensorField& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<edgeTensorField>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    return Laplacian;
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacian
(
    const tmp<edgeTensorField>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > Laplacian
    (
        fac::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}
*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
