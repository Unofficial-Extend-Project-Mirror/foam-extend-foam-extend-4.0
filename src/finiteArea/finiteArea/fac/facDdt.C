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

#include "facDdt.H"
#include "faMesh.H"
#include "faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ddt
(
    const dimensioned<Type> dt,
    const faMesh& mesh
)
{
    return fa::faDdtScheme<Type>::New
    (
        mesh,
        mesh.ddtScheme("ddt(" + dt.name() + ')')
    )().facDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ddt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + vf.name() + ')')
    )().facDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    )().facDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
ddt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    )().facDdt(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
