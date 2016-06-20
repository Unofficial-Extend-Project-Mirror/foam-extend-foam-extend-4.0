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

#include "fvcLaplacian.H"
#include "fvMesh.H"
#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, scalar>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().laplacianScheme(name)
    )().fvcLaplacian(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian(vf, "laplacian(" + vf.name() + ')');
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return gamma*fvc::laplacian(vf, name);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return gamma*fvc::laplacian
    (
        vf, "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().laplacianScheme(name)
    )().fvcLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian
    (
        tgamma,
        vf,
        "laplacian(" + tgamma().name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    return fvc::laplacian
    (
        gamma,
        tvf,
        "laplacian(" + gamma.name() + ',' + tvf().name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    return fvc::laplacian
    (
        tgamma,
        tvf,
        "laplacian(" + tgamma().name() + ',' + tvf().name() + ')'
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().laplacianScheme(name)
    )().fvcLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> > laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> > laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
