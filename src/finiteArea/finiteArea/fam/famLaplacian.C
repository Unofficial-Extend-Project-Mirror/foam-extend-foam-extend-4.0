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

#include "areaFields.H"
#include "edgeFields.H"
#include "faMatrix.H"
#include "faLaplacianScheme.H"
#include "edgeInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    edgeScalarField Gamma
    (
        IOobject
        (
            "gamma",
            vf.time().constant(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return fam::laplacian(Gamma, vf);
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    edgeScalarField Gamma
    (
        IOobject
        (
            "gamma",
            vf.time().constant(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return fam::laplacian(Gamma, vf, name);
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const dimensionedScalar& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    edgeScalarField Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fam::laplacian(Gamma, vf);
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const dimensionedScalar& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    edgeScalarField Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fam::laplacian(Gamma, vf, name);
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const areaScalarField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const areaScalarField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().famLaplacian(gamma, vf);
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const tmp<areaScalarField>& tgamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > Laplacian(fam::laplacian(tgamma(), vf));
    tgamma.clear();
    return Laplacian;
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const tmp<areaScalarField>& tgamma,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type> > Laplacian(fam::laplacian(tgamma(), vf, name));
    tgamma.clear();
    return Laplacian;
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const edgeScalarField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().famLaplacian(gamma, vf);
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const tmp<edgeScalarField>& tgamma,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type> > tLaplacian = fam::laplacian(tgamma(), vf, name);
    tgamma.clear();
    return tLaplacian;
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const edgeScalarField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}

template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const tmp<edgeScalarField>& tgamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam(fam::laplacian(tgamma(), vf));
    tgamma.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const areaTensorField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    return fam::laplacian
    (
        (mesh.Le() & fac::interpolate(gamma) & mesh.Le())
        /sqr(mesh.magLe()),
        vf
    );
}

template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const tmp<areaTensorField>& tgamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > Laplacian = fam::laplacian(tgamma(), vf);
    tgamma.clear();
    return Laplacian;
}


template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const edgeTensorField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    return fam::laplacian
    (
        (mesh.Le() & gamma & mesh.Le())/sqr(mesh.magLe()),
        vf
    );
}

template<class Type>
tmp<faMatrix<Type> >
laplacian
(
    const tmp<edgeTensorField>& tgamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > Laplacian = fam::laplacian(tgamma(), vf);
    tgamma.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
