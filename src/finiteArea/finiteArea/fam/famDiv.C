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

#include "famDiv.H"
#include "faMesh.H"
#include "faMatrix.H"
#include "faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> >
div
(
    const edgeScalarField& flux,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().divScheme(name)
    )().famDiv(flux, vf);
}

template<class Type>
tmp<faMatrix<Type> >
div
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type> > Div(fam::div(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<faMatrix<Type> >
div
(
    const edgeScalarField& flux,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::div(flux, vf, "div("+flux.name()+','+vf.name()+')');
}

template<class Type>
tmp<faMatrix<Type> >
div
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > Div(fam::div(tflux(), vf));
    tflux.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
