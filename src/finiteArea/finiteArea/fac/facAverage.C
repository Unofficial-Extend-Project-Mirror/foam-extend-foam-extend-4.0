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

#include "facAverage.H"
#include "facEdgeIntegrate.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
average
(
    const GeometricField<Type,  faePatchField, edgeMesh>& ssf
)
{
    const faMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, faPatchField, areaMesh> > taverage
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "average("+ssf.name()+')',
                ssf.instance(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            ssf.dimensions()
        )
    );

    GeometricField<Type, faPatchField, areaMesh>& av = taverage();

    av.internalField() = 
    (
        edgeSum(mesh.magLe()*ssf)/edgeSum(mesh.magLe())
    )().internalField();

    forAll(av.boundaryField(), patchi)
    {
        av.boundaryField()[patchi] = ssf.boundaryField()[patchi];
    }

    av.correctBoundaryConditions();

    return taverage;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
average
(
    const tmp<GeometricField<Type,  faePatchField, edgeMesh> >& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > taverage
    (
        fac::average(tssf())
    );
    tssf.clear();
    return taverage;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
average
(
    const GeometricField<Type, faPatchField, areaMesh>& vtf
)
{
    return fac::average(linearEdgeInterpolate(vtf));
}


template<class Type> 
tmp<GeometricField<Type, faPatchField, areaMesh> >
average
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvtf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > taverage
    (
        fac::average(tvtf())
    );
    tvtf.clear();
    return taverage;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
